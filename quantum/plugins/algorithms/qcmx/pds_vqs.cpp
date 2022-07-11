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
#include "pds_vqs.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/IndexedViewHelper.h"
#include "ObservableTransform.hpp"
#include "Optimizer.hpp"
#include "Utils.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "PauliOperator.hpp"
#include "xacc_observable.hpp"

#include <algorithm>
#include <complex>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <optional>

using namespace xacc;
namespace {
std::optional<xacc::quantum::PauliOperator>
getGenerator(const InstPtr &in_inst) {
  if (in_inst->name() == "Rx") {
    return xacc::quantum::PauliOperator({{in_inst->bits()[0], "X"}}, 0.5);
  }
  if (in_inst->name() == "Ry") {
    return xacc::quantum::PauliOperator({{in_inst->bits()[0], "Y"}}, 0.5);
  }
  if (in_inst->name() == "Rz") {
    return xacc::quantum::PauliOperator({{in_inst->bits()[0], "Z"}}, 0.5);
  }

  return std::nullopt;
}
} // namespace

using namespace xacc::quantum;

#ifdef _XACC_NO_STD_BETA
#include "boost/math/special_functions/beta.hpp"
namespace std {
double beta(double x, double y) { return boost::math::beta(x, y); }
} // namespace std
#endif

namespace xacc {
namespace algorithm {

bool PDS_VQS::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    xacc::error("Obs was false.");
    return false;
  }

  if (!parameters.pointerLikeExists<Accelerator>("accelerator")) {
    xacc::error("Acc was false.");
    return false;
  }

  if (!parameters.keyExists<int>("cmx-order")) {
    xacc::error("CMX order was false.");
    return false;
  }

  if (!parameters.pointerLikeExists<CompositeInstruction>("ansatz")) {
    xacc::error("Ansatz was false.");
    return false;
  }

  if (!parameters.pointerLikeExists<Optimizer>("optimizer")) {
    xacc::error("Optimizer was false.");
    return false;
  }

  accelerator = parameters.getPointerLike<Accelerator>("accelerator");
  order = parameters.get<int>("cmx-order");
  kernel = parameters.getPointerLike<CompositeInstruction>("ansatz");
  optimizer = parameters.getPointerLike<Optimizer>("optimizer");

  // check if Observable needs JW
  if (std::dynamic_pointer_cast<PauliOperator>(xacc::as_shared_ptr(
          parameters.getPointerLike<Observable>("observable")))) {
    observable = xacc::as_shared_ptr(
        parameters.getPointerLike<Observable>("observable"));
  } else {
    auto jw = xacc::getService<ObservableTransform>("jw");
    observable = jw->transform(xacc::as_shared_ptr(
        parameters.getPointerLike<Observable>("observable")));
  }

  if (parameters.keyExists<double>("threshold")) {
    threshold = parameters.get<double>("threshold");
    xacc::info("Ignoring measurements with coefficient below = " +
               std::to_string(threshold));
  }

  // metric
  if (parameters.stringExists("metric")) {
    metric = parameters.getString("metric");
  }

  // this is if we want to optimize excited states
  if (parameters.keyExists<int>("n-roots")) {
    nRoots = parameters.get<int>("n-roots");
  }

  return true;
}

const std::vector<std::string> PDS_VQS::requiredParameters() const {
  return {"observable", "accelerator", "cmx-order", "ansatz", "optimizer"};
}

void PDS_VQS::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  // First gather the operators for all required moments
  std::vector<std::shared_ptr<Observable>> momentOperators(2 * order - 1);
  auto momentOperator = PauliOperator("I");
  for (int i = 0; i < 2 * order - 1; i++) {
    momentOperator *= (*std::dynamic_pointer_cast<PauliOperator>(observable));
    momentOperators[i] = std::make_shared<PauliOperator>(momentOperator);
  }

  // now get the unique terms in all moments operators
  auto uniqueTerms = getUniqueTerms(momentOperators);

  std::vector<double> optSpectrum(order);
  double optEnergy = 10000; // just some ridiculous number

  OptFunction f(
      [&, this](const std::vector<double> &x, std::vector<double> &dx) {

        auto evaled = kernel->operator()(x);
        auto kernels = uniqueTerms->observe(evaled);
        auto nEnergyKernels = kernels.size();

        // now get circuits for gradients of unique terms
        auto parameterShift = xacc::getGradient("parameter-shift");
        parameterShift->initialize({{"observable", uniqueTerms}, {"shift-scalar", 0.5}});
        auto gradKernels = parameterShift->getGradientExecutions(
            xacc::as_shared_ptr(kernel), x);

        auto nGradKernels = gradKernels.size();
        kernels.insert(kernels.end(), gradKernels.begin(), gradKernels.end());

        // get circuits for metric
        auto nMetricKernels = 0;
        if (optimizer->name() == "riemann" && metric != "GD") {
          auto metricKernels =
              getMetricCircuits(xacc::as_shared_ptr(kernel), x);
          nMetricKernels = metricKernels.size();
          kernels.insert(kernels.end(), metricKernels.begin(),
                         metricKernels.end());
        }

        // excecute all circuits
        auto tmpBuffer = qalloc(buffer->size());
        accelerator->execute(tmpBuffer, kernels);

        // keep track of buffers
        auto buffers = tmpBuffer->getChildren();
        auto energyBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(
            buffers.begin(), buffers.begin() + nEnergyKernels);
        auto gradBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(
            buffers.begin() + nEnergyKernels,
            buffers.begin() + nEnergyKernels + nGradKernels);

        std::vector<std::shared_ptr<AcceleratorBuffer>> metricBuffers(
            nMetricKernels);
        if (optimizer->name() == "riemann" && metric != "GD") {
          metricBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(
              buffers.begin() + nEnergyKernels + nGradKernels, buffers.end());
        }

        // compute moments and spectrum
        auto moments = computeMoments(momentOperators, energyBuffers);
        auto spectrum = PDS(moments);

        // get ground state energy
        auto energy = *std::min_element(spectrum.begin(), spectrum.end());
        std::cout << "Energy: " << energy << "\n";
        if (energy < optEnergy) {
          optEnergy = energy;
          optSpectrum = spectrum;
        }

        // get denominator and derivative of the polynomial in eq. 17
        auto denominator = -order * std::pow(energy, order - 1);
        Eigen::VectorXd polynomialDerivative = Eigen::VectorXd::Zero(order);
        polynomialDerivative(order - 1) = 1.0;
        for (auto i = 0; i < order - 1; i++) {
          denominator -= (order - i - 1) * X(i) * std::pow(energy, order - i - 2);
          polynomialDerivative(i) = std::pow(energy, order - i - 1);
        }
        polynomialDerivative /= denominator;

        // get moments gradients
        auto nParams = x.size();
        auto nTerms = uniqueTerms->getNonIdentitySubTerms().size();
        auto momentsGradients = computeMomentsGradients(momentOperators,
                                                   gradBuffers, nTerms, nParams);

        // get the ground state energy gradient
        Eigen::VectorXd dE = Eigen::VectorXd::Zero(nParams);
        for (int i = 0; i < nParams; i++) {

          // get gradient wrt parameter
          auto gradMomentsParam = momentsGradients[i];

          // get dM
          Eigen::MatrixXd dM(order, order);
          dM = computeMatrixM(gradMomentsParam);
          // M(order, order) = <trial|trial> = 1
          // so it has no angle dependence, thus its derivative is 0
          dM(order - 1, order - 1) = 0.0;

          // get dY
          Eigen::VectorXd dY(order);
          dY = computeVectorY(gradMomentsParam);

          // get dX
          Eigen::VectorXd dX(order);
          dX = M.colPivHouseholderQr().solve(-dY - dM * X).transpose();

          // get dE
          dE(i) = polynomialDerivative.dot(dX);
        }

        if (optimizer->name() != "riemann") {
          dx = std::vector<double>(dE.data(), dE.data() + dE.size());
          std::cout << "x = " << x << " dx =" << dx << "\n";
          return optEnergy;
        }

        if (metric == "GD") {
          dx = std::vector<double>(dE.data(), dE.data() + dE.size());
        } else {
   
          Eigen::MatrixXd metricMatrix(nParams, nParams);
          metricMatrix = constructMetricTensorMatrix(metricBuffers, nParams);
          Eigen::VectorXd new_dx(nParams);
          new_dx = metricMatrix.colPivHouseholderQr().solve(dE).transpose();

          optimizer->appendOption(
              "new-dx", std::vector<double>(new_dx.data(),
                                             new_dx.data() + new_dx.size()));
          dx = std::vector<double>(new_dx.data(), new_dx.data() + new_dx.size());
        }

        return energy;
      },
      kernel->nVariables());

  auto result = optimizer->optimize(f);
  // add energy to the buffer
  buffer->addExtraInfo("spectrum", optSpectrum);
  buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
  buffer->addExtraInfo("opt-params", ExtraInfo(result.second));
  return;
}

// Compute energy from CMX
// J. Phys. A: Math. Gen. 17, 625 (1984)
// and Int. J. Mod. Phys. B 9, 2899 (1995)
std::vector<double> PDS_VQS::PDS(const std::vector<double> &moments) const {

  M = computeMatrixM(moments);
  Eigen::VectorXd Y(order);
  Y = computeVectorY(moments);

  // Solve linear system Ma = -b
  X = M.colPivHouseholderQr().solve(-Y).transpose();

  // Energy spectrum is given by roots of the polynomial
  // P_n(x) = \sum_0^n a_i x^(n - i)
  // Find roots by computing the companion matrix of P_n(x)
  // and diagonalizing it
  //
  // companionMatrix = |0    ... a_(n - 1)|
  //                   |I_n  ...    a_0   |
  //
  Eigen::MatrixXd companionMatrix = Eigen::MatrixXd::Zero(order, order);
  companionMatrix.bottomLeftCorner(order - 1, order - 1) =
      Eigen::MatrixXd::Identity(order - 1, order - 1);
  for (int i = 0; i < order; i++) {
    companionMatrix(i, order - 1) = -X(order - i - 1);
  }

  // get roots from eigenvalues
  // and return lowest energy eigenvalue
  Eigen::EigenSolver<Eigen::MatrixXd> es(companionMatrix);
  std::vector<double> spectrum(order);
  std::transform(es.eigenvalues().begin(), es.eigenvalues().end(),
                 spectrum.begin(),
                 [](std::complex<double> c) { return std::real(c); });

  return spectrum;
}

std::shared_ptr<Observable> PDS_VQS::getUniqueTerms(
    const std::vector<std::shared_ptr<Observable>> momentOps) const {

  std::vector<std::string> uniqueTermsStrings;
  auto uniqueTermsPtr = std::make_shared<PauliOperator>();
  for (auto &op : momentOps) {
    for (auto &subTerm : op->getNonIdentitySubTerms()) {
      auto termPtr = std::dynamic_pointer_cast<PauliOperator>(subTerm)->begin();
      auto nameTerm = termPtr->first;
      auto op = termPtr->second.ops();

      if (!xacc::container::contains(uniqueTermsStrings, nameTerm)) {
        uniqueTermsStrings.push_back(nameTerm);
        uniqueTermsPtr->operator+=(*std::make_shared<PauliOperator>(op));
      }
    }
  }

  return uniqueTermsPtr;
}

std::vector<double> PDS_VQS::computeMoments(
    const std::vector<std::shared_ptr<Observable>> momentOperators,
    const std::vector<std::shared_ptr<AcceleratorBuffer>> buffers) const {

  std::vector<double> moments;
  for (int i = 0; i < 2 * order - 1; i++) {
    double expval = 0.0;
    if (momentOperators[i]->getIdentitySubTerm()) {
      expval +=
          std::real(momentOperators[i]->getIdentitySubTerm()->coefficient());
    }

    for (auto subTerm : momentOperators[i]->getNonIdentitySubTerms()) {
      auto termName =
          std::dynamic_pointer_cast<PauliOperator>(subTerm)->begin()->first;

      for (auto buffer : buffers) {
        if (buffer->name() == termName) {
          expval += std::real(subTerm->coefficient() *
                              buffer->getExpectationValueZ());
          break;
        }
      }
    }
    moments.push_back(expval);
  }
  return moments;
}

Eigen::MatrixXd
PDS_VQS::computeMatrixM(const std::vector<double> &moments) const {

  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(order, order);
  for (int i = 0; i < order; i++) {
    for (int j = i; j < order; j++) {
      if (i + j == 2 * (order - 1)) {
        ret(i, j) = 1.0;
      } else {
        ret(i, j) = moments[2 * order - (i + j) - 3];
      }
      ret(j, i) = ret(i, j);
    }
  }

  return ret;
}

Eigen::VectorXd
PDS_VQS::computeVectorY(const std::vector<double> &moments) const {

  Eigen::VectorXd ret = Eigen::VectorXd::Zero(order);
  for (int i = 0; i < order; i++) {
    ret(i) = moments[2 * order - i - 2];
  }

  return ret;
}

std::vector<std::vector<double>> PDS_VQS::computeMomentsGradients(
    const std::vector<std::shared_ptr<Observable>> momentOperators,
    const std::vector<std::shared_ptr<AcceleratorBuffer>> buffers,
    const int nTerms, const int nParams) const {

  // loop over all parameters
  std::vector<std::vector<double>> momentsGradients;
  for (int paramIdx = 0; paramIdx < nParams; paramIdx++) {

    std::vector<std::pair<std::string, double>> paramUniqueTermGrad(nTerms);
    for (int termIdx = 0; termIdx < nTerms; termIdx++) {

      auto termName = buffers[termIdx + 2 * paramIdx * nTerms]->name();
      auto termGrad =
          (buffers[termIdx + 2 * paramIdx * nTerms]->getExpectationValueZ() -
           buffers[termIdx + (2 * paramIdx + 1) * nTerms]->getExpectationValueZ()) /
          2.0;
      paramUniqueTermGrad[termIdx] = {termName, termGrad};

    }

    std::vector<double> momentsGradientsParams(momentOperators.size());
    for (int m = 0; m < 2 * order -1; m++) {

      double momentGrad = 0.0;
      for (auto &subTerm : momentOperators[m]->getNonIdentitySubTerms()) {

        auto termName =
            std::dynamic_pointer_cast<PauliOperator>(subTerm)->begin()->first;
        auto coeff = subTerm->coefficient().real();

        for (auto &pair : paramUniqueTermGrad) {
          if (termName == pair.first) {
            momentGrad += coeff * pair.second;
            break;
          }
        }

      }

      momentsGradientsParams[m] = momentGrad;
    }

    momentsGradients.push_back(momentsGradientsParams);
  }
  return momentsGradients;
}

ObservedKernels
PDS_VQS::getMetricCircuits(std::shared_ptr<CompositeInstruction> in_circuit,
                           const std::vector<double> &in_x) const {

  // Layering the circuit:
  m_layers = ParametrizedCircuitLayer::toParametrizedLayers(in_circuit);
  auto m_nbParams = in_x.size();
  std::vector<std::string> paramNames;
  assert(in_circuit->getVariables().size() == in_x.size());
  for (const auto &param : in_circuit->getVariables()) {
    paramNames.emplace_back(param);
  }

  std::vector<std::shared_ptr<CompositeInstruction>> metricTensorKernels;
  for (auto &layer : m_layers) {
    auto kernels = constructMetricTensorSubCircuit(layer, paramNames, in_x);

    // Insert *non-identity* kernels only
    size_t kernelIdx = 0;
    const auto isIdentityTerm = [](xacc::quantum::PauliOperator &in_pauli) {
      return in_pauli.getNonIdentitySubTerms().size() == 0;
    };

    //if (metric == "NGD") {
      for (size_t i = 0; i < layer.kiTerms.size(); ++i) {
        if (!isIdentityTerm(layer.kiTerms[i])) {
          metricTensorKernels.emplace_back(kernels[kernelIdx]);
          m_metricTermToIdx.emplace(layer.kiTerms[i].toString(),
                                    metricTensorKernels.size() - 1);
        }
        kernelIdx++;
      //}
    }

    for (size_t i = 0; i < layer.kikjTerms.size(); ++i) {
      if (!isIdentityTerm(layer.kikjTerms[i])) {
        metricTensorKernels.emplace_back(kernels[kernelIdx]);
        m_metricTermToIdx.emplace(layer.kikjTerms[i].toString(),
                                  metricTensorKernels.size());
      }
      kernelIdx++;
    }
  }

  return metricTensorKernels;
}

ObservedKernels PDS_VQS::constructMetricTensorSubCircuit(
    ParametrizedCircuitLayer &io_layer,
    const std::vector<std::string> &in_varNames,
    const std::vector<double> &in_varVals) const {

  assert(in_varNames.size() == in_varVals.size());
  // We need to observe all the generators of parametrized gates plus products
  // of them.
  std::vector<xacc::quantum::PauliOperator> KiTerms;
  std::vector<xacc::quantum::PauliOperator> KiKjTerms;

  for (const auto &parOp : io_layer.ops) {
    auto opGenerator = getGenerator(parOp);
    assert(opGenerator.has_value());
    KiTerms.emplace_back(opGenerator.value());
  }

  for (const auto &genOp1 : KiTerms) {
    for (const auto &genOp2 : KiTerms) {
      auto kikjOp = genOp1 * genOp2;
      assert(kikjOp.nTerms() == 1);
      KiKjTerms.emplace_back(kikjOp);
    }
  }

  auto gateRegistry = xacc::getService<IRProvider>("quantum");
  auto circuitToObs = gateRegistry->createComposite("__LAYER__COMPOSITE__");
  std::vector<std::shared_ptr<xacc::CompositeInstruction>> obsComp;
  std::vector<double> resolvedParams;
  for (const auto &op : io_layer.preOps) {
    if (op->isParameterized() && op->getParameter(0).isVariable()) {
      const auto varName = op->getParameter(0).as<std::string>();
      circuitToObs->addVariable(varName);
      const auto iter =
          std::find(in_varNames.begin(), in_varNames.end(), varName);
      assert(iter != in_varNames.end());
      const size_t idx = std::distance(in_varNames.begin(), iter);
      assert(idx < in_varVals.size());
      resolvedParams.emplace_back(in_varVals[idx]);
    }
  }

  circuitToObs->addInstructions(io_layer.preOps);
  // Resolves the pre-ops now that we have passed that layer.
  auto resolvedCirc = resolvedParams.empty()
                          ? circuitToObs
                          : circuitToObs->operator()(resolvedParams);

  size_t NUM_KI_TERMS = 0;
  size_t NUM_KIKJ_TERMS = 0;
  //if (metric == "NGD") {
    io_layer.kiTerms = KiTerms;
    for (auto &term : KiTerms) {
      auto obsKernels = term.observe(resolvedCirc);
      assert(obsKernels.size() == 1);
      obsComp.emplace_back(obsKernels[0]);
    }
    NUM_KI_TERMS = io_layer.paramInds.size();
  //}

  io_layer.kikjTerms = KiKjTerms;
  for (auto &term : KiKjTerms) {
    auto obsKernels = term.observe(resolvedCirc);
    assert(obsKernels.size() == 1);
    obsComp.emplace_back(obsKernels[0]);
  }
  NUM_KIKJ_TERMS = io_layer.paramInds.size() * io_layer.paramInds.size();

  // Validate the expected count.
  assert(obsComp.size() == NUM_KI_TERMS + NUM_KIKJ_TERMS);
  return obsComp;
}

Eigen::MatrixXd PDS_VQS::constructMetricTensorMatrix(
    const std::vector<std::shared_ptr<xacc::AcceleratorBuffer>> &in_results,
    const int nParams) const {

  Eigen::MatrixXd gMat = Eigen::MatrixXd::Zero(nParams, nParams);
  size_t blockIdx = 0;
  const auto getExpectationForTerm =
      [&](xacc::quantum::PauliOperator &in_pauli) {
        if (m_metricTermToIdx.find(in_pauli.toString()) ==
            m_metricTermToIdx.end()) {
          return 1.0;
        }
        return in_results[m_metricTermToIdx[in_pauli.toString()]]
            ->getExpectationValueZ();
      };

  for (auto &layer : m_layers) {

    const auto nbParamsInBlock = layer.paramInds.size();
    // Constructs the block diagonal matrices
    Eigen::MatrixXd blockMat =
        Eigen::MatrixXd::Zero(nbParamsInBlock, nbParamsInBlock);

    for (size_t i = 0; i < nbParamsInBlock; ++i) {
      for (size_t j = 0; j < nbParamsInBlock; ++j) {

        // Entry = <KiKj> - <Ki><Kj>
        // second_order_ev[i, j] - first_order_ev[i] * first_order_ev[j]

        auto firstOrderTerm1 = layer.kiTerms[i];
        auto firstOrderTerm2 = layer.kiTerms[j];
        auto secondOrderTerm = layer.kiTerms[i] * layer.kiTerms[j];

        const double factor1 = firstOrderTerm1.coefficient().real();
        const double factor2 = firstOrderTerm2.coefficient().real();

        const double secondOrderTermExp =
            getExpectationForTerm(secondOrderTerm);
        double value = secondOrderTermExp;

        if (metric == "NGD") {
          const double firstOrderTerm1Exp =
              getExpectationForTerm(firstOrderTerm1);
          const double firstOrderTerm2Exp =
              getExpectationForTerm(firstOrderTerm2);
          value -= firstOrderTerm1Exp * firstOrderTerm2Exp;
        }

        // const double value = factor1*factor2*(secondOrderTermExp -
        // firstOrderTerm1Exp * firstOrderTerm2Exp);

        blockMat(i, j) = factor1 * factor2 * value;
      }
    }

    for (size_t i = 0; i < nbParamsInBlock; ++i) {
      for (size_t j = 0; j < nbParamsInBlock; ++j) {
        gMat(blockIdx + i, blockIdx + j) = blockMat(i, j);
      }
    }

    blockIdx += nbParamsInBlock;
  }

  return gMat;
}

std::vector<ParametrizedCircuitLayer>
ParametrizedCircuitLayer::toParametrizedLayers(
    const std::shared_ptr<xacc::CompositeInstruction> &in_circuit) {
  const auto variables = in_circuit->getVariables();
  if (variables.empty()) {
    xacc::error("Input circuit is not parametrized.");
  }

  // Assemble parametrized circuit layer:
  std::vector<ParametrizedCircuitLayer> layers;
  ParametrizedCircuitLayer currentLayer;
  std::set<size_t> qubitsInLayer;
  for (size_t instIdx = 0; instIdx < in_circuit->nInstructions(); ++instIdx) {
    const auto &inst = in_circuit->getInstruction(instIdx);
    if (!inst->isParameterized()) {
      if (currentLayer.ops.empty()) {
        currentLayer.preOps.emplace_back(inst);
      } else {
        currentLayer.postOps.emplace_back(inst);
      }
    } else {
      const auto &instParam = inst->getParameter(0);
      if (instParam.isVariable()) {
        const auto varName = instParam.as<std::string>();
        const auto iter =
            std::find(variables.begin(), variables.end(), varName);
        assert(iter != variables.end());
        assert(inst->bits().size() == 1);
        const auto bitIdx = inst->bits()[0];
        // This qubit line was already acted upon in this layer by a
        // parametrized gate. Hence, start a new layer.
        if (!currentLayer.postOps.empty() ||
            xacc::container::contains(qubitsInLayer, bitIdx)) {
          assert(!currentLayer.ops.empty() && !currentLayer.paramInds.empty());
          layers.emplace_back(std::move(currentLayer));
          qubitsInLayer.clear();
        }
        currentLayer.ops.emplace_back(inst);
        currentLayer.paramInds.emplace_back(
            std::distance(iter, variables.begin()));
        qubitsInLayer.emplace(bitIdx);
      } else {
        if (currentLayer.ops.empty()) {
          currentLayer.preOps.emplace_back(inst);
        } else {
          currentLayer.postOps.emplace_back(inst);
        }
      }
    }
  }

  assert(!currentLayer.ops.empty() && !currentLayer.paramInds.empty());
  layers.emplace_back(std::move(currentLayer));

  // Need to make a copy to do two rounds:
  auto layersCopied = layers;

  // Add pre-ops and post-ops:
  for (size_t i = 1; i < layers.size(); ++i) {
    // Pre-ops:
    auto &thisLayerPreOps = layers[i].preOps;
    const auto &previousLayerOps = layers[i - 1].ops;
    const auto &previousLayerPreOps = layers[i - 1].preOps;
    const auto &previousLayerPostOps = layers[i - 1].postOps;

    // Insert pre-ops, ops, and *raw* post-ops of the previous layer to the
    // begining
    thisLayerPreOps.insert(thisLayerPreOps.begin(),
                           previousLayerPostOps.begin(),
                           previousLayerPostOps.end());
    thisLayerPreOps.insert(thisLayerPreOps.begin(), previousLayerOps.begin(),
                           previousLayerOps.end());
    thisLayerPreOps.insert(thisLayerPreOps.begin(), previousLayerPreOps.begin(),
                           previousLayerPreOps.end());
  }

  for (int i = layers.size() - 2; i >= 0; --i) {
    // Post-ops:
    auto &thisLayerPostOps = layers[i].postOps;
    const auto &nextLayerOps = layers[i + 1].ops;
    const auto &nextLayerPostOps = layers[i + 1].postOps;
    // Get the original pre-ops
    const auto &nextLayerPreOps = layersCopied[i + 1].preOps;

    // Insert next layer pre-ops, ops and its post op
    thisLayerPostOps.insert(thisLayerPostOps.end(), nextLayerPreOps.begin(),
                            nextLayerPreOps.end());
    thisLayerPostOps.insert(thisLayerPostOps.end(), nextLayerOps.begin(),
                            nextLayerOps.end());
    thisLayerPostOps.insert(thisLayerPostOps.end(), nextLayerPostOps.begin(),
                            nextLayerPostOps.end());
  }

  // Verify:
  for (const auto &layer : layers) {
    assert(layer.preOps.size() + layer.ops.size() + layer.postOps.size() ==
           in_circuit->nInstructions());
    assert(!layer.paramInds.empty() && !layer.ops.empty());
  }

  return layers;
}

OptResult RiemannianMetricOptimizer::optimize(OptFunction &function) {
  auto dim = function.dimensions();
  double tol = 1e-6;
  int maxeval = 1000;

  std::vector<double> deltaX;
  if (options.keyExists<std::vector<double>>("delta-x")) {
    deltaX = options.get<std::vector<double>>("delta-x");
  }

  std::vector<double> x(dim, 0.0), dx(dim, 0.0);
  x[0] = xacc::constants::pi / 2.0;

  if (options.keyExists<std::vector<double>>("initial-parameters")) {
    x = options.get_with_throw<std::vector<double>>("initial-parameters");
  } else if (options.keyExists<std::vector<int>>("initial-parameters")) {
    auto tmpx = options.get<std::vector<int>>("initial-parameters");
    x = std::vector<double>(tmpx.begin(), tmpx.end());
  }

  if (options.keyExists<int>("maxeval")) {
    maxeval = options.get<int>("maxeval");
    xacc::info("max function evaluations set to " + std::to_string(maxeval));
  }

  double step = 0.05;
  if (options.keyExists<double>("step")) {
    maxeval = options.get<double>("step");
    xacc::info("step size set to " + std::to_string(maxeval));
  }

  int iter = 0;
  double error = 1.0, oldEnergy = 0.0;

  while (iter < maxeval || std::abs(error) > tol) {

    auto currentEnergy = function(x, dx);
    std::cout << dx << "  " << iter << "\n";
    error = currentEnergy - oldEnergy;
    oldEnergy = currentEnergy;

    std::transform(x.begin(), x.end(), deltaX.begin(), x.begin(),
                   std::plus<double>());

    std::cout << oldEnergy << "  " << iter << "\n";

    exit(0);
    iter++;
  }

  return OptResult{oldEnergy, x};
}

} // namespace algorithm
} // namespace xacc
