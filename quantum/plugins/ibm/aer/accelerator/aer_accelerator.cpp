/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Thien Nguyen - initial API and implementation
 *******************************************************************************/
#include "aer_accelerator.hpp"
#include "Accelerator.hpp"
#include "Utils.hpp"
#include "aer_noise_model.hpp"
#include "CommonGates.hpp"
#include "CountGatesOfTypeVisitor.hpp"
#include "InstructionIterator.hpp"

#include "controllers/aer_controller.hpp"
#include "controllers/controller_execute.hpp"
#include "noise/readout_error.hpp"

#include "cppmicroservices/BundleActivator.h"
#include "cppmicroservices/BundleContext.h"
#include "cppmicroservices/ServiceProperties.h"
#include "simulators/density_matrix/densitymatrix.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "aer_custatevec.hpp"

#include <bitset>
#include <optional>
#include "QObjGenerator.hpp"
#include "py-aer/aer_python_adapter.hpp"

namespace xacc {
namespace quantum {

#define IS_INTEGRAL(T)                                                         \
  typename std::enable_if<std::is_integral<T>::value>::type * = 0

template <class T>
std::string integral_to_binary_string(T byte, IS_INTEGRAL(T)) {
  std::bitset<sizeof(T) * 8> bs(byte);
  return bs.to_string();
}
std::string hex_string_to_binary_string(const std::string &hex) {
  auto hex_char_to_bin = [](char c) -> std::string {
    switch (toupper(c)) {
    case '0':
      return "0000";
    case '1':
      return "0001";
    case '2':
      return "0010";
    case '3':
      return "0011";
    case '4':
      return "0100";
    case '5':
      return "0101";
    case '6':
      return "0110";
    case '7':
      return "0111";
    case '8':
      return "1000";
    case '9':
      return "1001";
    case 'A':
      return "1010";
    case 'B':
      return "1011";
    case 'C':
      return "1100";
    case 'D':
      return "1101";
    case 'E':
      return "1110";
    case 'F':
      return "1111";
    }
    throw std::runtime_error("Invalid Hex character!");
    return "";
  };

  // start with '0x'
  assert(hex.substr(0, 2) == "0x");
  std::string bin;
  for (int i = 2; i != hex.length(); ++i) {
    bin += hex_char_to_bin(hex[i]);
  }

  return bin;
}

HeterogeneousMap AerAccelerator::getProperties() {
  auto props = physical_backend_properties;
  // Insert 'shots' data
  if (m_simtype == "qasm") {
    props.insert("shots", m_shots);
  }
  return props;
}

void AerAccelerator::initialize(const HeterogeneousMap &params) {

  m_options = params;
  noise_model.clear();
  m_simtype = "qasm";
  connectivity.clear();

  void* state = getState();

  xacc_to_qobj = xacc::getCompiler("qobj");
  if (params.keyExists<int>("shots")) {
    m_shots = params.get<int>("shots");
  }

  if (params.keyExists<int>("seed")) {
    m_seed = params.get<int>("seed");
  }
  if (params.stringExists("sim-type")) {
    if (!xacc::container::contains(
            std::vector<std::string>{"qasm", "statevector", "pulse",
                                     "density_matrix", "matrix_product_state"},
            params.getString("sim-type"))) {
      xacc::warning("[Aer] warning, invalid sim-type (" +
                    params.getString("sim-type") +
                    "), must be qasm or statevector or density_matrix or "
                    "matrix_product_state.");
    } else {
      m_simtype = params.getString("sim-type");
    }
  }

  // A custom sim method is specified with 'shots';
  // use qasm sim but set the method in QObj to force the method.
  // e.g., statevector and density_matrix will not do shots, just returns the
  // 'state'.
  if (m_simtype != "qasm" && params.keyExists<int>("shots")) {
    m_simtype = "qasm";
  }

  if (params.stringExists("backend")) {
    auto ibm_noise_model = xacc::getService<NoiseModel>("IBM");
    ibm_noise_model->initialize(params);
    auto json_str = ibm_noise_model->toJson();
    noise_model = nlohmann::json::parse(json_str);
    auto ibm = xacc::getAccelerator("ibm:" + params.getString("backend"));
    physical_backend_properties = ibm->getProperties();
    if (m_simtype == "pulse") {
      // If pulse mode, must contribute pulse cmd-def.
      ibm->contributeInstructions();
    }
    // Set connectivity based on the backend:
    connectivity = ibm->getConnectivity();
  } else if (params.stringExists("noise-model")) {
    std::string noise_model_str = params.getString("noise-model");
    // Check if this is a file name
    std::ifstream test(noise_model_str);
    if (test) {
      std::string str((std::istreambuf_iterator<char>(test)),
                      std::istreambuf_iterator<char>());
      noise_model = nlohmann::json::parse(str);
    } else {
      noise_model = nlohmann::json::parse(params.getString("noise-model"));
    }
    
    // need this for ro error decorator
    std::vector<double> p01s, p10s;
    auto errors = noise_model["errors"];
    for (auto error : errors) {
      if (error["operations"].get<std::vector<std::string>>() ==
          std::vector<std::string>{"measure"}) {
        auto probs =
            error["probabilities"].get<std::vector<std::vector<double>>>();
        auto meas1prep0 = probs[0][1];
        auto meas0prep1 = probs[1][0];
        p01s.push_back(meas0prep1);
        p10s.push_back(meas1prep0);
      }
    }
    physical_backend_properties.insert("p01s", p01s);
    physical_backend_properties.insert("p10s", p10s);

  }
  AER::Hacks::maybe_load_openmp("");
  initialized = true;
}
double AerAccelerator::calcExpectationValueZ(
    const std::vector<std::pair<double, double>> &in_stateVec,
    const std::vector<std::size_t> &in_bits) {
  const auto hasEvenParity =
      [](size_t x, const std::vector<size_t> &in_qubitIndices) -> bool {
    size_t count = 0;
    for (const auto &bitIdx : in_qubitIndices) {
      if (x & (1ULL << bitIdx)) {
        count++;
      }
    }
    return (count % 2) == 0;
  };

  double result = 0.0;
  for (uint64_t i = 0; i < in_stateVec.size(); ++i) {
    result += (hasEvenParity(i, in_bits) ? 1.0 : -1.0) *
              std::norm(std::complex<double>(in_stateVec[i].first,
                                             in_stateVec[i].second));
  }

  return result;
}

double AerAccelerator::calcExpectationValueZFromDensityMatrix(
    const std::vector<std::vector<std::pair<double, double>>> &in_densityMat,
    const std::vector<std::size_t> &in_bits) {
  const auto hasEvenParity =
      [](size_t x, const std::vector<size_t> &in_qubitIndices) -> bool {
    size_t count = 0;
    for (const auto &bitIdx : in_qubitIndices) {
      if (x & (1ULL << bitIdx)) {
        count++;
      }
    }
    return (count % 2) == 0;
  };

  double result = 0.0;
  for (uint64_t i = 0; i < in_densityMat.size(); ++i) {
    const auto &diag_elem = in_densityMat[i][i];
    // The diag. elements of the DM should be real.
    assert(std::abs(diag_elem.second) < 1e-3);
    // When using DM elements, don't need the square-norm.
    result += ((hasEvenParity(i, in_bits) ? 1.0 : -1.0) * diag_elem.first);
  }

  return result;
}

void AerAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> program) {
  if (m_simtype == "qasm") {
    auto qobj_str = xacc_to_qobj->translate(program);

    nlohmann::json j = nlohmann::json::parse(qobj_str)["qObject"];
    j["config"]["shots"] = m_shots;
    j["config"]["noise_model"] = noise_model;
    // If a seed was set:
    if (m_seed > 0) {
      j["config"]["seed_simulator"] = m_seed;
    }

    if (m_options.stringExists("sim-type")) {
      const std::string requestedMethod = m_options.getString("sim-type");
      if (requestedMethod == "statevector") {
        j["config"]["method"] = "statevector";
      } else if (requestedMethod == "density_matrix") {
        j["config"]["method"] = "density_matrix";
      } else if (requestedMethod == "matrix_product_state") {
        j["config"]["method"] = "matrix_product_state";
      } else {
        j["config"]["method"] = "automatic";
      }
    }

    // xacc::set_verbose(true);
    // xacc::info("Shots Qobj:\n" + j.dump(2));
    // Initialize QOBJ
    AER::Qobj qobj(j);
    AER::Controller controller;
    // Set config
    controller.set_config(qobj.config);
    assert(qobj.circuits.size() == 1);

    if (m_options.keyExists<std::vector<std::complex<double>>>(
            "initial_state")) {
      const std::vector<std::complex<double>> intial_state =
          m_options.get<std::vector<std::complex<double>>>("initial_state");
      if (intial_state.size() != (1ULL << buffer->size())) {
        xacc::error("Dimension mismatch in 'initial_state', expected size of " +
                    std::to_string(1ULL << buffer->size()));
      }
      const double norm = std::accumulate(
          intial_state.begin(), intial_state.end(), 0.0,
          [](double current_norm, const std::complex<double> &val) {
            return current_norm + std::norm(val);
          });
      if (std::abs(norm - 1.0) > 1e-12) {
        xacc::error("Initial state vector input is not normalized.");
      }
      // Add initialize op to set the initial input state
      AER::Operations::Op op;
      op.type = AER::Operations::OpType::initialize;
      op.name = "initialize";
      for (int q = 0; q < buffer->size(); ++q) {
        op.qubits.emplace_back(q);
      }
      op.params = intial_state;
      qobj.circuits[0].ops.insert(qobj.circuits[0].ops.begin(), op);
      qobj.circuits[0].set_params(false);
    }

    // Run qobj circuits
    auto result = controller.execute(qobj.circuits, qobj.noise_model, qobj.config);
    if (result.status != AER::Result::Status::completed) {
      xacc::error("Failed to complete the simulation! Error: " +
                  result.message);
    }
    assert(result.results.size() == 1);
    auto results = result.results.front();
    auto counts =
        results.data.to_json()["counts"].get<std::map<std::string, int>>();
    // Process bitStr to be an n-Measure string in msb
    CountGatesOfTypeVisitor<Measure> cc(program);
    const int nMeasures = cc.countGates();
    for (auto &kv : counts) {
      std::string hexStr = kv.first;
      int nOccurrences = kv.second;
      auto bitStr = hex_string_to_binary_string(hexStr);
      std::string actual(nMeasures, '0');
      for (int i = 0; i < nMeasures; i++) {
        if (bitStr.length() > i) {
          actual[actual.length() - 1 - i] = bitStr[bitStr.length() - i - 1];
        }
      }
      buffer->appendMeasurement(actual, nOccurrences);
    }
  } else if (m_simtype == "pulse") {
    // Get the correct QObject Generator
    auto qobjGen = xacc::getService<QObjGenerator>("pulse");
    auto chosenBackend = nlohmann::json::parse(physical_backend_properties.getString("config-json"));
    // Unused
    const std::string getBackendPropsResponse = "{}";
    auto defaults_response = nlohmann::json::parse(physical_backend_properties.getString("defaults-json"));
    auto ibmPulseAssembler = xacc::getService<IRTransformation>("ibm-pulse");
    auto kernel = xacc::ir::asComposite(program->clone());
    
    // Remove measures, add measure all (at the end)
    kernel->clear();
    InstructionIterator iter(program);
    std::vector<size_t> measured_bits;
    while (iter.hasNext()) {
      auto next = iter.next();
      if (!next->isComposite() && next->name() != "Measure") {
        kernel->addInstruction(next);
      } else if (next->name() == "Measure") {
        measured_bits.push_back(next->bits()[0]);
      }
    }

    auto provider = xacc::getIRProvider("quantum");
    for (size_t qId = 0; qId < buffer->size(); ++qId) {
      kernel->addInstruction(provider->createInstruction("Measure", qId));
    }

    // Assemble pulse composite from the input kernel.
    ibmPulseAssembler->apply(kernel, nullptr);    
    // Generate the QObject JSON
    auto qobjJsonStr = qobjGen->getQObjJsonStr({kernel}, m_shots, chosenBackend,
                                           getBackendPropsResponse,
                                           connectivity, defaults_response);
    xacc::info("Qobj:\n" + qobjJsonStr);
    auto hamiltonianJson = chosenBackend["hamiltonian"];
    // Remove unrelated fields (could contain problematic characters)
    hamiltonianJson.erase("description");
    hamiltonianJson.erase("h_latex"); 
    xacc::info("Hamiltonian Json:\n" + hamiltonianJson.dump());
    const auto dt = chosenBackend["dt"].get<double>();
    const auto qubitFreqEst = defaults_response["qubit_freq_est"].get<std::vector<double>>();
    const auto uLoFreqs = chosenBackend["u_channel_lo"];
    std::vector<int> uLoRefs;
    for (auto loIter = uLoFreqs.begin(); loIter != uLoFreqs.end(); ++loIter) {
      auto uLoConfig = *((*loIter).begin());
      uLoRefs.emplace_back(uLoConfig["q"].get<int>());
    }
    // Run the simulation via Python
    const std::string resultJson =
        xacc::aer::runPulseSim(hamiltonianJson.dump(), dt, qubitFreqEst, uLoRefs, qobjJsonStr);
    auto result_json = nlohmann::json::parse(resultJson);
    auto count_json = result_json["counts"].get<std::map<std::string, int>>();
    for (const auto &[hexStr, nOccurrences] : count_json) {
      auto bitStr = hex_string_to_binary_string(hexStr);
      // Process bitStr to be an n-Measure string in msb
      std::string actual = "";
      const auto nMeasures = measured_bits.size();
      for (int i = 0; i < nMeasures; i++) {
        actual += "0";
      }

      for (int i = 0; i < nMeasures; i++) {
        actual[actual.length() - 1 - i] = bitStr[bitStr.length() - measured_bits[i] - 1];
      }

      buffer->appendMeasurement(actual, nOccurrences);
    }

    auto state_vector = result_json["statevector"].get<std::vector<std::pair<double, double>>>();
    buffer->addExtraInfo("state", state_vector);
  } else if (m_simtype == "density_matrix") {
    // remove all measures, don't need them
    auto tmp = xacc::ir::asComposite(program->clone());
    tmp->clear();
    AER::reg_t measured_bits;
    InstructionIterator iter(program);
    while (iter.hasNext()) {
      auto next = iter.next();
      if (!next->isComposite() && next->name() != "Measure") {
        tmp->addInstruction(next);
      } else if (next->name() == "Measure") {
        measured_bits.push_back(next->bits()[0]);
      }
    }
    // In "density_matrix" simulation mode, always include Id gates
    // if they are explicitly added to the Composite.
    auto qobj_str = xacc_to_qobj->translate(tmp, {{"skip-id-gates", false}});    
    nlohmann::json qObjJson = nlohmann::json::parse(qobj_str)["qObject"];
    qObjJson["config"]["method"] = "density_matrix";
    qObjJson["config"]["enable_truncation"] = false;
    qObjJson["config"]["noise_model"] = noise_model;
    AER::DensityMatrix::State<AER::QV::DensityMatrix<double>> densityMat;
    AER::RngEngine rng;
    AER::Qobj qobj(qObjJson);
    // std::cout << "QObj:\n" << qobj_str << "\n";
    assert(qobj.circuits.size() == 1);
    auto circ = qobj.circuits[0];
    // Output data container
    AER::ExperimentResult data;
    densityMat.initialize_creg(circ.num_memory, circ.num_registers);
    densityMat.initialize_qreg(buffer->size());
    if (m_options.keyExists<std::vector<std::complex<double>>>(
            "initial_state")) {
      const std::vector<std::complex<double>> intial_state =
          m_options.get<std::vector<std::complex<double>>>("initial_state");
      if (intial_state.size() != (1ULL << buffer->size())) {
        xacc::error("Dimension mismatch in 'initial_state', expected size of " +
                    std::to_string(1ULL << buffer->size()));
      }
      const double norm = std::accumulate(
          intial_state.begin(), intial_state.end(), 0.0,
          [](double current_norm, const std::complex<double> &val) {
            return current_norm + std::norm(val);
          });
      if (std::abs(norm - 1.0) > 1e-12) {
        xacc::error("Initial state vector input is not normalized.");
      }
      densityMat.qreg().initialize_from_vector(intial_state);
    }
    // std::cout << "Num op: " << circ.ops.size() << "\n";
    if (!noise_model.empty()) {
      auto noise = qobj.noise_model;
      noise.enable_superop_method();
      auto opt_circ = noise.sample_noise(
          circ, rng, AER::Noise::NoiseModel::Method::superop);
      densityMat.apply_ops(opt_circ.ops.begin(), opt_circ.ops.end(), data, rng);
    } else {
      densityMat.apply_ops(circ.ops.begin(), circ.ops.end(), data, rng);
    }
    // std::cout << "Result: \n" << data.to_json().dump() << "\n";
    const double exp_val = densityMat.qreg().expval_pauli(
        measured_bits, std::string(measured_bits.size(), 'Z'));

    auto dmData = densityMat.move_to_matrix(0);
    const auto length = dmData.size();
    assert(length == (1ULL << buffer->size()) * (1ULL << buffer->size()));
    std::complex<double> *dmPtr = dmData.move_to_buffer();
    std::vector<std::pair<double, double>> flattenDm;
    flattenDm.reserve(length);
    for (int i = 0; i < length; ++i) {
      flattenDm.emplace_back(dmPtr[i].real(), dmPtr[i].imag());
    }
    buffer->addExtraInfo("density_matrix", flattenDm);

    ExecutionInfo::DensityMatrixType dm(1ULL << buffer->size());
    auto *startAddr = dmPtr;
    for (auto &row : dm) {
      auto *endAddr = startAddr + (1ULL << buffer->size());
      row.assign(startAddr, endAddr);
      startAddr = endAddr;
    }
    m_executionInfo = {
        {ExecutionInfo::DmKey,
         std::make_shared<ExecutionInfo::DensityMatrixType>(std::move(dm))}};
    buffer->addExtraInfo("exp-val-z", exp_val);
  }
  else {
    assert(m_simtype == "statevector");
    // statevector
    // remove all measures, don't need them
    auto tmp = xacc::ir::asComposite(program->clone());
    tmp->clear();
    AER::reg_t measured_bits;
    InstructionIterator iter(program);
    while (iter.hasNext()) {
      auto next = iter.next();
      if (!next->isComposite() && next->name() != "Measure") {
        tmp->addInstruction(next);
      } else if (next->name() == "Measure") {
        measured_bits.push_back(next->bits()[0]);
      }
    }
    auto qobj_str = xacc_to_qobj->translate(tmp);
    nlohmann::json qObjJson = nlohmann::json::parse(qobj_str)["qObject"];
    qObjJson["config"]["method"] = "statevector";
    qObjJson["config"]["enable_truncation"] = false;
    AER::Statevector::State<AER::QV::QubitVector<double>> stateVec;
    AER::RngEngine rng;
    AER::Qobj qobj(qObjJson);
    // std::cout << "QObj:\n" << qobj_str << "\n";
    assert(qobj.circuits.size() == 1);
    auto circ = qobj.circuits[0];
    // Output data container
    AER::ExperimentResult data;
    stateVec.initialize_creg(circ.num_memory, circ.num_registers);
    stateVec.initialize_qreg(buffer->size());
    if (m_options.keyExists<std::vector<std::complex<double>>>(
            "initial_state")) {
      const std::vector<std::complex<double>> intial_state =
          m_options.get<std::vector<std::complex<double>>>("initial_state");
      if (intial_state.size() != (1ULL << buffer->size())) {
        xacc::error("Dimension mismatch in 'initial_state', expected size of " +
                    std::to_string(1ULL << buffer->size()));
      }
      const double norm = std::accumulate(
          intial_state.begin(), intial_state.end(), 0.0,
          [](double current_norm, const std::complex<double> &val) {
            return current_norm + std::norm(val);
          });
      if (std::abs(norm - 1.0) > 1e-12) {
        xacc::error("Initial state vector input is not normalized.");
      }
      stateVec.qreg().initialize_from_vector(intial_state);
    }
    // std::cout << "Num op: " << circ.ops.size() << "\n";
    stateVec.apply_ops(circ.ops.begin(), circ.ops.end(), data, rng);
    // std::cout << "Result: \n" << data.to_json().dump() << "\n";
    const double exp_val = stateVec.qreg().expval_pauli(
        measured_bits, std::string(measured_bits.size(), 'Z'));

    auto stateVecData = stateVec.move_to_vector(0);
    const auto length = stateVecData.size();
    std::complex<double> *statePtr = stateVecData.move_to_buffer();
    ExecutionInfo::WaveFuncType wavefn(statePtr, statePtr + length);
    // std::cout << "HOWDY: " << wavefn.size() << "\n";
    // for (int i = 0; i < wavefn.size(); ++i) {
    //   std::cout << wavefn[i] << "\n";
    // }
    buffer->addExtraInfo("exp-val-z", exp_val);
    m_executionInfo = {
        {ExecutionInfo::WaveFuncKey,
         std::make_shared<ExecutionInfo::WaveFuncType>(std::move(wavefn))}};
  }
}

void AerAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>>
        compositeInstructions) {
  for (auto &f : compositeInstructions) {
    auto tmpBuffer =
        std::make_shared<xacc::AcceleratorBuffer>(f->name(), buffer->size());
    execute(tmpBuffer, f);
    buffer->appendChild(f->name(), tmpBuffer);
  }
}

void AerAccelerator::apply(std::shared_ptr<AcceleratorBuffer> buffer,
                           std::shared_ptr<Instruction> inst) {
  static AER::Statevector::State<AER::QV::QubitVector<double>> stateVec;
  static auto provider = xacc::getIRProvider("quantum");
  static AER::RngEngine rng;
  static std::set<AcceleratorBuffer*> knownBuffers;
  if (!xacc::container::contains(knownBuffers, buffer.get())) {
    stateVec.initialize_qreg(buffer->size());
    knownBuffers.emplace(buffer.get());
  }
  if (inst->isComposite() || inst->isAnalog()) {
    xacc::error("Only gates are allowed.");
  }
  xacc::info("Apply: " + inst->toString());
  if (inst->name() != "Measure") {
    auto tempComp = provider->createComposite("tmp");
    tempComp->addInstruction(inst);
    auto qobj_str = xacc_to_qobj->translate(tempComp);
    auto qObjJson = nlohmann::json::parse(qobj_str)["qObject"];
    qObjJson["config"]["enable_truncation"] = false;
    AER::Qobj qobj(qObjJson);
    // std::cout << "QObj:\n" << qobj_str << "\n";
    assert(qobj.circuits.size() == 1);
    auto circ = qobj.circuits[0];
    // Output data container
    AER::ExperimentResult data;
    // std::cout << "Num op: " << circ.ops.size() << "\n";
    stateVec.apply_ops(circ.ops.begin(), circ.ops.end(), data, rng);
    // std::cout << "Result: \n" << data.to_json().dump() << "\n";
    // auto wavefn = stateVec.copy_to_vector(0);
    // std::cout << "HOWDY: " << wavefn.size() << "\n";
    // for (int i = 0; i < wavefn.size(); ++i) {
    //   std::cout << wavefn[i] << "\n";
    // }
  }
  //   stateVec.add_creg_to_data(data);
  // If it was a Measure op:
  else {
    AER::reg_t measured_bits{inst->bits()[0]};
    const auto probs = stateVec.qreg().probabilities(measured_bits);
    assert(probs.size() == 2);
    assert(std::abs(1.0 - probs[0] - probs[1]) < 1e-12);
    // Randomly pick outcome and return pair
    auto outcome = rng.rand_int(probs);
    // std::cout << "Outcome: " << outcome << "\n";
    assert(outcome == 0 || outcome == 1);
    // std::cout << "Probs: " << probs[0] << "  " << probs[1] << "\n";
    std::vector<std::complex<double>> mdiag(2, 0.);
    mdiag[outcome] = 1. / std::sqrt(probs[outcome]);
    stateVec.qreg().apply_diagonal_matrix(measured_bits, mdiag);
    buffer->measure(inst->bits()[0], outcome);
  }
}

void IbmqNoiseModel::initialize(const HeterogeneousMap &params) {
  m_nbQubits = 0;
  m_qubitT1.clear();
  m_qubitT2.clear();
  m_gateErrors.clear();
  m_gateDurations.clear();
  m_roErrors.clear();
  m_connectivity.clear();
  // Note: we support both remote backend JSON and cache JSON string.
  // So that we can test this with offline JSON.
  if (params.stringExists("backend") || params.stringExists("backend-json")) {
    const auto backendJsonStr = [&]() {
      if (params.stringExists("backend")) {
        auto ibm = xacc::getAccelerator("ibm:" + params.getString("backend"));
        auto props = ibm->getProperties().get<std::string>("total-json");
        m_connectivity = ibm->getConnectivity();
        return props;
      }
      return params.getString("backend-json");
    }();

    auto backEndJson = nlohmann::json::parse(backendJsonStr);
    // Cache the backend properties JSON.
    m_backendPropertiesJson = backendJsonStr;
    // Parse qubit data:
    auto qubitsData = backEndJson["qubits"];
    size_t nbQubit = 0;
    for (auto qubitIter = qubitsData.begin(); qubitIter != qubitsData.end();
         ++qubitIter) {
      std::optional<double> meas0Prep1, meas1Prep0;
      // Each qubit contains an array of properties.
      for (auto probIter = qubitIter->begin(); probIter != qubitIter->end();
           ++probIter) {
        const auto probObj = *probIter;
        const std::string probName = probObj["name"].get<std::string>();
        const double probVal = probObj["value"].get<double>();
        const std::string unit = probObj["unit"].get<std::string>();
        if (probName == "T1") {
          assert(unit == "µs" || unit == "us" || unit == "ns");
          m_qubitT1.emplace_back(
              (unit == "µs" || unit == "us") ? 1000.0 * probVal : probVal);
        }

        if (probName == "T2") {
          assert(unit == "µs" || unit == "us" || unit == "ns");
          m_qubitT2.emplace_back(
              (unit == "µs" || unit == "us") ? 1000.0 * probVal : probVal);
        }

        if (probName == "prob_meas0_prep1") {
          assert(unit.empty());
          meas0Prep1 = probVal;
        }

        if (probName == "prob_meas1_prep0") {
          assert(unit.empty());
          meas1Prep0 = probVal;
        }
      }
      assert(meas0Prep1.has_value() && meas1Prep0.has_value());
      m_roErrors.emplace_back(std::make_pair(*meas0Prep1, *meas1Prep0));

      nbQubit++;
    }
    m_nbQubits = nbQubit;
    assert(m_qubitT1.size() == m_nbQubits);
    assert(m_qubitT2.size() == m_nbQubits);
    assert(m_roErrors.size() == m_nbQubits);

    // Parse gate data:
    auto gateData = backEndJson["gates"];
    for (auto gateIter = gateData.begin(); gateIter != gateData.end();
         ++gateIter) {
      auto gateObj = *gateIter;
      const std::string gateName = gateObj["name"].get<std::string>();
      auto gateParams = gateObj["parameters"];
      for (auto it = gateParams.begin(); it != gateParams.end(); ++it) {
        auto paramObj = *it;
        const std::string paramName = paramObj["name"].get<std::string>();
        if (paramName == "gate_length") {
          const std::string unit = paramObj["unit"].get<std::string>();
          assert(unit == "µs" || unit == "us" || unit == "ns");
          const double gateLength =
              (unit == "µs" || unit == "us")
                  ? 1000.0 * paramObj["value"].get<double>()
                  : paramObj["value"].get<double>();
          const bool insertOk =
              m_gateDurations.insert({gateName, gateLength}).second;
          if (gateName.rfind("sx", 0) == 0) {
            // This device properties were specified in the { rz, sx, cx } basis.
            // We add equivalent data for { u1, u2, u3 } basis set as well.
            // For every sx gate, add equiv. data for u1, u2, u3 gates.
            // Note: some of our noisy simulators don't split a gate into multiple sx gates,
            // hence, they cannot simulate noise at that level.
            const std::string qubitOperandSuffix = gateName.substr(2);
            const std::string u1GateName = "u1" + qubitOperandSuffix;
            const std::string u2GateName = "u2" + qubitOperandSuffix;
            const std::string u3GateName = "u3" + qubitOperandSuffix;
            m_gateDurations.insert({u1GateName, 0.0});
            m_gateDurations.insert({u2GateName, gateLength});
            m_gateDurations.insert({u3GateName, 2.0 * gateLength});
          }
          // Must not contain duplicates.
          assert(insertOk);
        }

        if (paramName == "gate_error") {
          assert(paramObj["unit"].get<std::string>().empty());
          const bool insertOk =
              m_gateErrors.insert({gateName, paramObj["value"].get<double>()})
                  .second;
          // Must not contain duplicates.
          assert(insertOk);
          if (gateName.rfind("sx", 0) == 0) {
            const std::string qubitOperandSuffix = gateName.substr(2);
            const std::string u1GateName = "u1" + qubitOperandSuffix;
            const std::string u2GateName = "u2" + qubitOperandSuffix;
            const std::string u3GateName = "u3" + qubitOperandSuffix;
            m_gateErrors.insert({u1GateName, 0.0});
            m_gateErrors.insert({u2GateName, paramObj["value"].get<double>()});
            m_gateErrors.insert({u3GateName, 2.0 * paramObj["value"].get<double>()});
          }
        }
      }
    }
  }
}

std::string
IbmqNoiseModel::getUniversalGateEquiv(xacc::quantum::Gate &in_gate) const {
  if (in_gate.bits().size() == 1 && in_gate.name() != "Measure") {
    // Note: rotation around Z is a noiseless *u1* operation;
    // *u2* operations are those that requires a half-length rotation;
    // *u3* operations are those that requires a full-length rotation.
    static const std::unordered_map<std::string, std::string>
        SINGLE_QUBIT_GATE_MAP{{"X", "u3"},   {"Y", "u3"},  {"Z", "u1"},
                              {"H", "u2"},   {"U", "u3"},  {"T", "u1"},
                              {"Tdg", "u1"}, {"S", "u1"},  {"Sdg", "u1"},
                              {"Rz", "u1"},  {"Rx", "u3"}, {"Ry", "u3"}};
    const auto iter = SINGLE_QUBIT_GATE_MAP.find(in_gate.name());
    // If cannot find the gate, just treat that as a noiseless u1 op.
    const std::string universalGateName =
        (iter == SINGLE_QUBIT_GATE_MAP.end()) ? "u1" : iter->second;
    return universalGateName + std::to_string(in_gate.bits()[0]);
  }

  if (in_gate.bits().size() == 2) {
    return "cx" + std::to_string(in_gate.bits()[0]) + "_" +
           std::to_string(in_gate.bits()[1]);
  }

  return "id" + std::to_string(in_gate.bits()[0]);
}

// Return gate time, T1, T2
std::tuple<double, double, double>
IbmqNoiseModel::relaxationParams(xacc::quantum::Gate &in_gate,
                                 size_t in_qubitIdx) const {
  const std::string universalGateName = getUniversalGateEquiv(in_gate);
  const auto gateDurationIter = m_gateDurations.find(universalGateName);
  assert(gateDurationIter != m_gateDurations.end());
  const double gateDuration = gateDurationIter->second;
  const double qubitT1 = m_qubitT1[in_qubitIdx];
  const double qubitT2 = m_qubitT2[in_qubitIdx];
  return std::make_tuple(gateDuration, qubitT1, qubitT2);
}

std::vector<std::vector<std::complex<double>>>
IbmqNoiseModel::thermalRelaxationChoiMat(double in_gateTime, double in_T1,
                                         double in_T2) const {
  const double rate1 = 1.0 / in_T1;
  const double pReset = 1.0 - std::exp(-in_gateTime * rate1);

  const double rate2 = 1.0 / in_T2;
  const double expT2 = std::exp(-in_gateTime * rate2);
  const double p0 = 1.0;
  const double p1 = 0.0;
  return {{1 - p1 * pReset, 0, 0, expT2},
          {0, p1 * pReset, 0, 0},
          {0, 0, p0 * pReset, 0},
          {expT2, 0, 0, 1 - p0 * pReset}};
}

std::vector<double> IbmqNoiseModel::calculateDepolarizing(
    xacc::quantum::Gate &in_gate,
    const std::vector<std::vector<std::complex<double>>> &in_relaxationError)
    const {
  //  Compute the depolarizing channel error parameter in the
  //  presence of T1/T2 thermal relaxation.
  //  Hence we have that the depolarizing error probability
  //  for the composed depolarization channel is
  //  p = dim * (F(E_relax) - F) / (dim * F(E_relax) - 1)
  const double averageThermalError =
      1.0 - averageGateFidelity(in_gate, in_relaxationError);
  const std::string universalGateName = getUniversalGateEquiv(in_gate);
  // Retrieve the error rate:
  const auto gateErrorIter = m_gateErrors.find(universalGateName);
  const double gateError =
      (gateErrorIter == m_gateErrors.end()) ? 0.0 : gateErrorIter->second;
  // If the backend gate error (estimated by randomized benchmarking) is more
  // than thermal relaxation error. We need to add depolarization to simulate
  // the total gate error.
  if (gateError > averageThermalError) {
    // Model gate error entirely as depolarizing error
    const double depolError = 2 * (gateError - averageThermalError) /
                              (2 * (1 - averageThermalError) - 1);
    return {depolError};
  }
  return {0.0};
}

double IbmqNoiseModel::averageGateFidelity(
    xacc::quantum::Gate &in_gate,
    const std::vector<std::vector<std::complex<double>>> &in_relaxationError)
    const {
  // We only handle single-qubit gates for now.
  if (in_gate.bits().size() == 1 && in_relaxationError.size() == 4) {
    // Note: this is based on the simplified amplitude-damping channel
    const double processFidelity =
        (in_relaxationError[0][0].real() + in_relaxationError[0][3].real() +
         in_relaxationError[3][0].real() + in_relaxationError[3][3].real()) /
        4.0;

    assert(processFidelity <= 1.0);
    const double averageFidelity = (4.0 * processFidelity + 1.0) / 5.0;
    return averageFidelity;
  }
  return 1.0;
}

double IbmqNoiseModel::gateErrorProb(xacc::quantum::Gate &gate) const {
  const auto universalGateName = getUniversalGateEquiv(gate);
  const auto gateErrorIter = m_gateErrors.find(universalGateName);
  return (gateErrorIter == m_gateErrors.end()) ? 0.0 : gateErrorIter->second;
}

std::vector<NoiseChannelKraus>
IbmqNoiseModel::getNoiseChannels(xacc::quantum::Gate &gate) const {
  std::vector<NoiseChannelKraus> krausOps;
  const auto noiseUtils = xacc::getService<NoiseModelUtils>("default");
  if (gate.bits().size() == 1 && gate.name() != "Measure") {
    // Amplitude damping + dephasing
    const auto [gateDuration, qubitT1, qubitT2] =
        relaxationParams(gate, gate.bits()[0]);
    const auto relaxationError =
        thermalRelaxationChoiMat(gateDuration, qubitT1, qubitT2);

    // Depolarization
    const auto dpAmpl = calculateDepolarizing(gate, relaxationError);
    if (!dpAmpl.empty()) {
      const double probDP = dpAmpl[0];
      const std::vector<std::vector<std::complex<double>>> depolErrorChoi{
          {1.0 - probDP / 2.0, 0., 0., 1.0 - probDP},
          {0., probDP / 2.0, 0., 0.},
          {0., 0., probDP / 2.0, 0.},
          {1.0 - probDP, 0., 0., 1.0 - probDP / 2.0}};
      const auto noiseUtils = xacc::getService<NoiseModelUtils>("default");
      krausOps.emplace_back(NoiseChannelKraus(gate.bits(), noiseUtils->choiToKraus(depolErrorChoi), KrausMatBitOrder::MSB));
    }
  }
  // For two-qubit gates, we currently only support
  // amplitude damping on both qubits (scaled by gate time).
  // We don't have ability to handle multi-qubit depolarization
  // in TNQVM yet.
  if (gate.bits().size() == 2) {

    for (const auto &qubitIdx : gate.bits()) {
      const auto [gateDuration, qubitT1, qubitT2] =
          relaxationParams(gate, qubitIdx);
      const auto relaxationError =
          thermalRelaxationChoiMat(gateDuration, qubitT1, qubitT2);
      krausOps.emplace_back(NoiseChannelKraus({qubitIdx}, noiseUtils->choiToKraus(relaxationError), KrausMatBitOrder::MSB));
    }
  }
  return krausOps;
}

std::string IbmqNoiseModel::toJson() const {
  // First, check if Qiskit is available,
  // try to use Qiskit util if possible.
  const std::string noiseModelJson =
      xacc::aer::noiseModelFromBackendProperties(m_backendPropertiesJson);
  if (!noiseModelJson.empty()) {
    xacc::info("Qiskit generated noise model:\n" + noiseModelJson);
    return noiseModelJson;
  }
  
  // No Qiskit, create the noise model by ourselves from the backend properties.
  // Aer noise model Json
  nlohmann::json noiseModel;
  std::vector<nlohmann::json> noiseElements;
  // Adds RO errors:
  const auto roErrors = readoutErrors();
  for (size_t qIdx = 0; qIdx < roErrors.size(); ++qIdx) {
    const auto &[meas0Prep1, meas1Prep0] = roErrors[qIdx];
    const std::vector<std::vector<double>> probs{{1 - meas1Prep0, meas1Prep0},
                                                 {meas0Prep1, 1 - meas0Prep1}};
    nlohmann::json element;
    element["type"] = "roerror";
    element["operations"] = std::vector<std::string>{"measure"};
    element["probabilities"] = probs;
    element["gate_qubits"] = std::vector<std::vector<std::size_t>>{{qIdx}};
    noiseElements.push_back(element);
  }

  const auto noiseUtils = xacc::getService<NoiseModelUtils>("default");
  // Add Kraus noise:
  // (1) Single-qubit gate noise:
  // Note: we must add noise ops for u2, u3, and cx gates:
  for (size_t qIdx = 0; qIdx < roErrors.size(); ++qIdx) {
    // For mapping purposes:
    // U2 == Hadamard gate
    // U3 == X gate
    Hadamard gateU2(qIdx);
    X gateU3(qIdx);
    const std::unordered_map<std::string, xacc::quantum::Gate *> gateMap{
        {"u2", &gateU2}, {"u3", &gateU3}};

    for (const auto &[gateName, gate] : gateMap) {
      const auto errorChannels = getNoiseChannels(*gate);
      nlohmann::json element;
      element["type"] = "qerror";
      element["operations"] = std::vector<std::string>{gateName};
      element["gate_qubits"] = std::vector<std::vector<std::size_t>>{{qIdx}};
      std::vector<nlohmann::json> krausOps;
      for (const auto &error : errorChannels) {
        const auto krausOpMats = error.mats;
        nlohmann::json instruction;
        instruction["name"] = "kraus";
        instruction["qubits"] = std::vector<std::size_t>{0};
        instruction["params"] = krausOpMats;
        krausOps.emplace_back(instruction);
      }
      element["instructions"] =
          std::vector<std::vector<nlohmann::json>>{krausOps};
      element["probabilities"] = std::vector<double>{1.0};
      noiseElements.push_back(element);
    }
  }

  // (2) Two-qubit gate noise:
  for (const auto &[qubit1, qubit2] : m_connectivity) {
    // We need to add noise for both CNOT directions
    // Note: the duration of them can be different,
    // hence the noise channels.
    CNOT cx1(std::vector<std::size_t>{(size_t)qubit1, (size_t)qubit2});
    CNOT cx2(std::vector<std::size_t>{(size_t)qubit2, (size_t)qubit1});
    const std::vector<xacc::quantum::Gate *> cxGates{&cx1, &cx2};
    for (const auto &cx : cxGates) {
      const auto errorChannels = getNoiseChannels(*cx);
      assert(errorChannels.size() == 2);
      nlohmann::json element;
      element["type"] = "qerror";
      element["operations"] =
          std::vector<std::string>{getUniversalGateEquiv(*cx)};
      element["gate_qubits"] =
          std::vector<std::vector<std::size_t>>{cx->bits()};
      std::vector<nlohmann::json> krausOps;
      size_t noiseBitIdx = 0;
      for (const auto &error : errorChannels) {
        const auto krausOpMats = error.mats;
        nlohmann::json instruction;
        instruction["name"] = "kraus";
        instruction["qubits"] = std::vector<std::size_t>{noiseBitIdx++};
        instruction["params"] = krausOpMats;
        krausOps.emplace_back(instruction);
      }
      assert(krausOps.size() == 2);
      element["instructions"] =
          std::vector<std::vector<nlohmann::json>>{krausOps};
      element["probabilities"] = std::vector<double>{1.0};
      noiseElements.push_back(element);
    }
  }

  noiseModel["errors"] = noiseElements;
  return noiseModel.dump(6);
}

std::vector<double> IbmqNoiseModel::averageSingleQubitGateFidelity() const {
  std::vector<double> result;
  // Use U3 gate fidelity:
  for (size_t qId = 0; qId < m_nbQubits; ++qId) {
    const std::string gateName = "u3_" + std::to_string(qId);
    const auto gateErrorIter = m_gateErrors.find(gateName);
    const double gateError =
        (gateErrorIter == m_gateErrors.end()) ? 0.0 : gateErrorIter->second;
    result.emplace_back(1.0 - gateError);
  }
  return result;
}

std::vector<std::tuple<size_t, size_t, double>>
IbmqNoiseModel::averageTwoQubitGateFidelity() const {
  std::vector<std::tuple<size_t, size_t, double>> result;
  for (const auto &[gateName, gateError] : m_gateErrors) {
    if (gateName.rfind("cx", 0) == 0) {
      // CNOT gate:
      const std::size_t pos = gateName.find("_");
      const std::string firstArg = gateName.substr(2, pos - 2);
      const std::string secondArg = gateName.substr(pos + 1);
      const auto firstBit = std::atoi(firstArg.c_str());
      const auto secondBit = std::atoi(secondArg.c_str());
      result.emplace_back(
          std::make_tuple(firstBit, secondBit, 1.0 - gateError));
    }
  }

  return result;
}

std::string
AerAccelerator::getNativeCode(std::shared_ptr<CompositeInstruction> program,
                              const HeterogeneousMap &config) {
  auto qobj_str = xacc_to_qobj->translate(program);
  nlohmann::json j = nlohmann::json::parse(qobj_str)["qObject"];
  return j.dump(2);
}
} // namespace quantum
} // namespace xacc

using namespace cppmicroservices;

namespace {

/**
 */
class US_ABI_LOCAL AerActivator : public BundleActivator {

public:
  AerActivator() {}

  /**
   */
  void Start(BundleContext context) {
    auto xt = std::make_shared<xacc::quantum::AerAccelerator>();
    context.RegisterService<xacc::Accelerator>(xt);
    context.RegisterService<xacc::NoiseModel>(
        std::make_shared<xacc::quantum::IbmqNoiseModel>());
  }

  /**
   */
  void Stop(BundleContext /*context*/) {}
};

} // namespace

CPPMICROSERVICES_EXPORT_BUNDLE_ACTIVATOR(AerActivator)
