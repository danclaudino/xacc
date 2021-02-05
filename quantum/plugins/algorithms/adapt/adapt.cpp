/*******************************************************************************
 * Copyright (c) 2020 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Daniel Claudino - initial API and implementation
 ******************************************************************************/

#include "adapt.hpp"
#include "FermionOperator.hpp"
#include "ObservableTransform.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include "Circuit.hpp"
#include "AlgorithmGradientStrategy.hpp"

#include <complex>
#include <memory>
#include <iomanip>
#include <sstream>
#include <string>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {
namespace algorithm {

bool ADAPT::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    xacc::info("Obs was false\n");
    return false;
  }

  if (!parameters.pointerLikeExists<Accelerator>("accelerator")) {
    xacc::info("Acc was false\n");
    return false;
  }

  if (!parameters.stringExists("pool")) {
    xacc::info("Pool was false\n");
    return false;
  }

  if (!parameters.stringExists("sub-algorithm")) {
    xacc::info("Sub algorithm was false\n");
    return false;
  }

  optimizer = parameters.get<std::shared_ptr<Optimizer>>("optimizer");
  observable = parameters.get<std::shared_ptr<Observable>>("observable");
  accelerator = parameters.get<std::shared_ptr<Accelerator>>("accelerator");
  pool = xacc::getService<OperatorPool>(parameters.getString("pool"));
  subAlgo = parameters.getString("sub-algorithm");

  if (parameters.keyExists<int>("maxiter")) {
    _maxIter = parameters.get<int>("maxiter");
  }

  if (parameters.keyExists<double>("adapt-threshold")) {
    _adaptThreshold = parameters.get<double>("adapt-threshold");
  }

  if (parameters.keyExists<double>("print-threshold")) {
    _printThreshold = parameters.get<double>("print-threshold");
  }

  if (parameters.keyExists<bool>("print-operators")) {
    _printOps = parameters.get<bool>("print-operators");
  }

  if (parameters.stringExists("gradient_strategy")) {
    gradStrategyName = parameters.getString("gradient_strategy");
  }

  if (parameters.pointerLikeExists<CompositeInstruction>("initial-state")) {
    initialState =
        parameters.get<std::shared_ptr<CompositeInstruction>>("initial-state");
  }

  if (parameters.keyExists<int>("n-electrons")) {
    _nElectrons = parameters.get<int>("n-electrons");
  }

  if ((subAlgo == "vqe" && !initialState) &&
      (subAlgo == "vqe" && !parameters.keyExists<int>("n-electrons"))) {

    xacc::info("VQE requires number of electrons or initial state.");
  }

  if (parameters.getString("pool") == "singlet-adapted-uccsd" &&
      parameters.keyExists<int>("n-electrons")) {

    pool->optionalParameters({std::make_pair("n-electrons", _nElectrons)});
  }

  // Check if Observable is Fermion or Pauli and manipulate accordingly
  //
  // if string has ^, it's FermionOperator
  if (observable->toString().find("^") != std::string::npos) {

    auto jw = xacc::getService<ObservableTransform>("jw");
    if (std::dynamic_pointer_cast<FermionOperator>(observable)) {
      observable = jw->transform(observable);
    } else {
      auto fermionObservable =
          xacc::quantum::getObservable("fermion", observable->toString());
      observable = jw->transform(
          std::dynamic_pointer_cast<Observable>(fermionObservable));
    }

    // observable is PauliOperator, but does not cast down to it
    // Not sure about the likelihood of this happening, but want to cover all
    // bases
  } else if (observable->toString().find("X") != std::string::npos ||
             observable->toString().find("Y") != std::string::npos ||
             observable->toString().find("Z") != std::string::npos &&
                 !std::dynamic_pointer_cast<PauliOperator>(observable)) {

    auto pauliObservable =
        xacc::quantum::getObservable("pauli", observable->toString());
    observable = std::dynamic_pointer_cast<Observable>(pauliObservable);
  }

  if (parameters.keyExists<std::vector<double>>("initial-parameters")) {
    initialParams = parameters.get<std::vector<double>>("initial-parameters");
  }

  if (parameters.keyExists<std::vector<int>>("initial-operators")) {
    initialOps = parameters.get<std::vector<int>>("initial-operators");
  }

  if (parameters.keyExists<double>("measurement-threshold")) {
    _measurementThreshold = parameters.get<double>("measurement-threshold");
    xacc::info("Ignoring measurements with coefficient below = " + std::to_string(_measurementThreshold));
  }

  _parameters = parameters;

  return true;
}

const std::vector<std::string> ADAPT::requiredParameters() const {
  return {"observable", "optimizer", "accelerator", "pool", "sub-algorithm"};
}

void ADAPT::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  // stream to direct printing
  std::stringstream ss;
  // ansatz parameters
  std::vector<double> x;
  // ansatz operator indices according to the pool numbering
  std::vector<int> ansatzOps;

  // Instantiate adaptive circuit and add instructions from initial state
  // initial state for chemistry VQE is usually HF
  // initial state for QAOA is all |+>
  auto provider = xacc::getIRProvider("quantum");
  auto ansatz = provider->createComposite("ansatz");

  // Generate operators in the pool
  auto operators = pool->generate(buffer->size());

  // construct initial state
  if (initialState) {

    for (auto &inst : initialState->getInstructions()) {
      ansatz->addInstruction(inst);
    }

  } else if (!initialParams.empty() && !initialOps.empty()) {

    xacc::info("Starting from provided ansatz operators and parameters");
    x = initialParams;
    ansatzOps = initialOps;

    for (int i = 0; i < ansatzOps.size(); i++) {
      ansatz->addVariable("x" + std::to_string(i));
      for (auto &inst :
           pool->getOperatorInstructions(ansatzOps[i], i)->getInstructions()) {
        ansatz->addInstruction(inst);
      }
    }

  } else {

    if (subAlgo == "vqe") {
      // Define the initial state, usually HF for chemistry problems
      std::size_t j;
      for (int i = 0; i < _nElectrons / 2; i++) {
        j = (std::size_t)i;
        ansatz->addInstruction(provider->createInstruction("X", {j}));
        j = (std::size_t)(i + buffer->size() / 2);
        ansatz->addInstruction(provider->createInstruction("X", {j}));
      }
    }

    if (subAlgo == "QAOA") {
      std::size_t j;
      for (int i = 0; i < buffer->size(); i++) {
        j = (std::size_t)i;
        ansatz->addInstruction(provider->createInstruction("H", {j}));
      }
    }
  }

  // Vector of commutators, need to compute them only once
  std::vector<std::shared_ptr<Observable>> commutators;
  if (subAlgo == "QAOA") {

    std::complex<double> i(0, 1);
    for (auto op : operators) {
      auto comm = observable->commutator(op);
      auto &tmp = *std::dynamic_pointer_cast<PauliOperator>(comm);
      tmp = tmp * i;
      commutators.push_back(std::make_shared<PauliOperator>(tmp));
    }

  } else {

    for (auto op : operators) {
      commutators.push_back(observable->commutator(op));
    }
  }

  xacc::info("Operator pool: " + pool->name());
  xacc::info("Number of operators in the pool: " +
             std::to_string(operators.size()));

  // energy from the previous iteration
  double oldEnergy = 0.0;

  // number of circuit implementations
  // lower bound to the number of measurements
  // # measurements = nImplCircuits * shots;
  int nImplCircuits = 0;

  // start ADAPT loop
  for (int iter = x.size(); iter < _maxIter; iter++) {

    xacc::info("Iteration: " + std::to_string(iter + 1));
    xacc::info("Computing [H, A]");
    ss << "Printing commutators with absolute value above "  << _printThreshold;
    xacc::info(ss.str());
    ss.str(std::string());

    if (subAlgo == "QAOA") {

      auto costHamiltonianGates = std::dynamic_pointer_cast<quantum::Circuit>(
          xacc::getService<Instruction>("exp_i_theta"));

      // Create instruction for new operator
      costHamiltonianGates->expand(
          {{"pauli", observable->toString()},
           {"param_id", "x" + std::to_string(ansatz->nVariables())}});

      ansatz->addVariable("x" + std::to_string(ansatz->nVariables()));
      for (auto &inst : costHamiltonianGates->getInstructions()) {
        ansatz->addInstruction(inst);
      }
      x.insert(x.begin(), 0.01);
    }

    int maxCommutatorIdx = 0;
    double maxCommutator = 0.0, gradientNorm = 0.0;

    // Loop over non-vanishing commutators
    // select the one with largest magnitude
    for (int opIdx = 0; opIdx < commutators.size(); opIdx++) {

      // only compute commutators if they aren't zero
      int nTermsCommutator =
          std::dynamic_pointer_cast<PauliOperator>(commutators[opIdx])
              ->getTerms()
              .size();
      if (nTermsCommutator != 0) {
        // observe the commutators with the updated circuit ansatz
        auto commutatorValue =
            measureOperator(commutators[opIdx], ansatz, x, buffer->size());

        if (abs(commutatorValue) > _printThreshold) {
          // Print number of instructions for computing <observable>
          xacc::info("Number of instructions for commutator calculation: " +
                     std::to_string(nTermsCommutator));
          ss << std::setprecision(12) << "[H," << opIdx
             << "] = " << commutatorValue;
          xacc::info(ss.str());
          ss.str(std::string());
        }

        // update maxCommutator
        if (abs(commutatorValue) > abs(maxCommutator)) {
          maxCommutatorIdx = opIdx;
          maxCommutator = commutatorValue;
        }

        // update gradient norm
        gradientNorm += commutatorValue * commutatorValue;

        // increment nImplCircuits
        nImplCircuits += nTermsCommutator;
      }

      // clear cachedMeasurements
      cachedMeasurements.clear();
    }

    ss << std::setprecision(12) << "Max gradient component: [H, "
       << maxCommutatorIdx << "] = " << maxCommutator << " a.u.";
    xacc::info(ss.str());
    ss.str(std::string());

    gradientNorm = std::sqrt(gradientNorm);
    ss << std::setprecision(12) << "Norm of gradient vector: " << gradientNorm
       << " a.u.";
    xacc::info(ss.str());
    ss.str(std::string());

    if (gradientNorm < _adaptThreshold) { // ADAPT converged

      xacc::info("ADAPT-" + subAlgo + " converged in " + std::to_string(iter) +
                 " iterations.");
      return;

    } else if (iter < _maxIter) { // Add operator and reoptimize

      xacc::info(subAlgo + " optimization of current ansatz.");

      // keep track of growing ansatz
      ansatzOps.push_back(maxCommutatorIdx);

      // Instruction service for the operator to be added to the ansatz
      auto maxCommutatorGate =
          pool->getOperatorInstructions(maxCommutatorIdx, ansatz->nVariables());

      // Label for new variable and add it to the circuit
      ansatz->addVariable("x" + std::to_string(ansatz->nVariables()));

      // Append new instructions to current circuit
      for (auto &inst : maxCommutatorGate->getInstructions()) {
        ansatz->addInstruction(inst);
      }

      // update vector of parameters
      double param;
      if (subAlgo == "vqe") {
        param = newParameter(operators[maxCommutatorIdx], ansatz, x, oldEnergy);
        x.insert(x.begin(), param);
        nImplCircuits += 2 * observable->getNonIdentitySubTerms().size();
      }
      if (subAlgo == "QAOA") {
        x.insert(x.begin(), xacc::constants::pi / 2.0);
      }
      optimizer->appendOption("initial-parameters", x);

      // Instantiate gradient class
      std::shared_ptr<AlgorithmGradientStrategy> gradientStrategy;
      if (!gradStrategyName.empty() && optimizer->isGradientBased()) {
        gradientStrategy =
            xacc::getService<AlgorithmGradientStrategy>(gradStrategyName);
        gradientStrategy->initialize({{"observable", observable}});
      }

      // Start subAlgo optimization
      _parameters.insert("ansatz", ansatz);
      _parameters.insert("observable", observable);
      auto sub_opt = xacc::getAlgorithm(subAlgo, _parameters);

      // if this is the first iteration, param above is already optimal
      double newEnergy;
      if (iter == 0) {
        newEnergy = sub_opt->execute(buffer, {param})[0];
      } else {
        sub_opt->execute(buffer);
        newEnergy = (*buffer)["opt-val"].as<double>();
        x = (*buffer)["opt-params"].as<std::vector<double>>();
      }
      nImplCircuits += observable->getNonIdentitySubTerms().size();
      oldEnergy = newEnergy;

      ss << std::setprecision(12) << "Energy at ADAPT iteration " << iter + 1
         << ": " << newEnergy;
      xacc::info(ss.str());
      ss.str(std::string());

      ss << std::setprecision(12) << "Parameters at ADAPT iteration "
         << iter + 1 << ": ";
      for (auto param : x) {
        ss << param << " ";
      }
      xacc::info(ss.str());
      ss.str(std::string());

      ss << "Ansatz at ADAPT iteration " << std::to_string(iter + 1) << ": ";
      for (auto op : ansatzOps) {
        ss << op << " ";
      }
      xacc::info(ss.str());
      ss.str(std::string());

      ss << "Circuit depth at ADAPT iteration " << std::to_string(iter + 1)
         << ": " << std::to_string(ansatz->depth());
      xacc::info(ss.str());
      ss.str(std::string());

      ss << "Total number of gates at ADAPT iteration "
         << std::to_string(iter + 1) << ": "
         << std::to_string(ansatz->nInstructions());
      xacc::info(ss.str());
      ss.str(std::string());

      ss << "Total number of implemented circuits at ADAPT iteration "
         << std::to_string(iter + 1) << ": " << std::to_string(nImplCircuits);
      xacc::info(ss.str());
      ss.str(std::string());

      // persisting XASM circuit string to the buffer
      auto xasm = xacc::getCompiler("xasm");
      auto xasmCode = xasm->translate(ansatz->operator()(x));

      // persisted relevant info to the buffer
      buffer->addExtraInfo("opt-val", ExtraInfo(oldEnergy));
      buffer->addExtraInfo("opt-params", ExtraInfo(x));
      buffer->addExtraInfo("opt-ansatz", ExtraInfo(ansatzOps));
      buffer->addExtraInfo("n-implemented-circuits", ExtraInfo(nImplCircuits));
      buffer->addExtraInfo("opt-circuit", ExtraInfo(xasmCode));

    } else {
      xacc::info("ADAPT-" + subAlgo + " did not converge in " +
                 std::to_string(_maxIter) + " iterations.");
      return;
    }
  }
}

// initializes new parameter from partial tomography
// where all previous parameters are frozen
// much better than starting from 0
double ADAPT::newParameter(const std::shared_ptr<Observable> op,
                           const std::shared_ptr<CompositeInstruction> kernel,
                           const std::vector<double> &x,
                           double zeroExpValue) const {

  // First we need to get the coefficient for the operator
  auto terms = op->getSubTerms();
  auto term = std::dynamic_pointer_cast<PauliOperator>(terms[0])
                  ->getTerms()
                  .begin()
                  ->second;
  double coeff;
  if (std::fabs(std::get<0>(term)) == 1.0) {
    coeff = 0.5;
  } else {
    coeff = 2.0 * std::fabs(std::get<0>(term));
  }

  // now we compute the expectation value for different angles
  // and get the tomography coefficients
  auto vqe = xacc::getAlgorithm("vqe", {{"observable", observable},
                                        {"accelerator", accelerator},
                                        {"ansatz", kernel}});

  std::vector<double> tmp_x;
  // if this is the first iteration, we don't have the
  // the energy/expectation value E(0.0, ...)
  // for the other iterations, we cache the oldEnergy
  if (zeroExpValue == 0.0) {
    tmp_x = x;
    tmp_x.insert(tmp_x.begin(), 0.0);
    zeroExpValue = vqe->execute(xacc::qalloc(observable->nBits()), tmp_x)[0];
  }

  tmp_x = x;
  tmp_x.insert(tmp_x.begin(), coeff * xacc::constants::pi / 2.0);
  auto plusExpValue = vqe->execute(xacc::qalloc(observable->nBits()), tmp_x)[0];

  tmp_x = x;
  tmp_x.insert(tmp_x.begin(), -coeff * xacc::constants::pi / 2.0);
  auto minusExpValue =
      vqe->execute(xacc::qalloc(observable->nBits()), tmp_x)[0];

  auto B = zeroExpValue - (plusExpValue + minusExpValue) / 2.0;
  auto C = (plusExpValue - minusExpValue) / 2.0;

  return std::atan2(-C, -B) * coeff;
}

double
ADAPT::measureOperator(const std::shared_ptr<Observable> obs,
                       const std::shared_ptr<CompositeInstruction> kernel,
                       const std::vector<double> &x,
                       const int bufferSize) const {

  // observe
  auto evaled = kernel->operator()(x);
  auto kernels = obs->observe(evaled);

  // we loop over all measured circuits
  // and check if that term has been measured
  // if so, we just multiply the measurement by the coefficient
  // We gather all new circuits into fsToExec and execute
  // Because these are all commutators, we don't need to worry about the I term
  std::vector<std::shared_ptr<CompositeInstruction>> fsToExec;
  std::vector<std::complex<double>> coefficients;
  double total = 0.0;
  for (auto &f : kernels) {
    std::complex<double> coeff = f->getCoefficient();
    if (cachedMeasurements.find(f->name()) != cachedMeasurements.end()) {
      total += std::real(coeff * cachedMeasurements[f->name()]);
    } else if (fabs(coeff) >= _measurementThreshold) {
      fsToExec.push_back(f);
      coefficients.push_back(coeff);
    }
  }

  // for circuits that have not been executed previously
  auto tmpBuffer = xacc::qalloc(bufferSize);
  accelerator->execute(tmpBuffer, fsToExec);
  auto buffers = tmpBuffer->getChildren();
  for (int i = 0; i < fsToExec.size(); i++) {
    auto expval = buffers[i]->getExpectationValueZ();
    total += std::real(expval * coefficients[i]);
    cachedMeasurements.emplace(fsToExec[i]->name(), expval);
  }

  return total;
}

} // namespace algorithm
} // namespace xacc