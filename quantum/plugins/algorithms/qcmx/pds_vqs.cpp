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
#include "ObservableTransform.hpp"
#include "Optimizer.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include "OperatorPool.hpp"

#include <algorithm>
#include <complex>
#include <iomanip>
#include <memory>
#include <sstream>

using namespace xacc;
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
  initialState = parameters.getPointerLike<CompositeInstruction>("ansatz");
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
    printThreshold = parameters.get<double>("threshold");
    xacc::info("Ignoring measurements with coefficient below = " +
               std::to_string(printThreshold));
  }

  if (parameters.stringExists("pool")) {
    poolName = parameters.getString("pool");
  }

  if (parameters.keyExists<int>("n-electrons")) {
    nElectrons = parameters.get<int>("n-electrons");
  }

  if (parameters.keyExists<bool>("adapt")) {
    adapt = parameters.get<bool>("adapt");
  }

  if (parameters.keyExists<int>("maxiter")) {
    maxIterations = parameters.get<int>("maxiter");
  }

  if (parameters.keyExists<double>("adapt-threshold")) {
    _adaptThreshold = parameters.get<double>("adapt-threshold");
  }

  return true;
}

const std::vector<std::string> PDS_VQS::requiredParameters() const {
  return {"observable", "accelerator", "cmx-order", "ansatz", "optimizer"};
}

void PDS_VQS::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  std::stringstream ss;
  ss << std::setprecision(12);

  // First gather the operators for all required moments
  std::vector<std::shared_ptr<Observable>> momentOperators(2 * order - 1);
  auto momentOperator = PauliOperator("I");
  for (int i = 0; i < 2 * order - 1; i++) {
    momentOperator *= (*std::dynamic_pointer_cast<PauliOperator>(observable));
    momentOperators[i] = std::make_shared<PauliOperator>(momentOperator);
  }

  auto provider = xacc::getIRProvider("quantum");
  auto kernel = provider->createComposite("ansatz");
  if (initialState->nVariables()) {
    kernel->addVariables(initialState->getVariables());
  }
  for (auto &inst : initialState->getInstructions()) {
    kernel->addInstruction(inst);
  }

  // now get the unique terms in all moments operators
  auto uniqueTerms = getUniqueTerms(momentOperators);

  std::vector<double> optSpectrum(order);
  double optEnergy = 10000; // just some ridiculous number

  // performance metrics
  int nTotalCircuits = 0, nTotalMeasurements = 0;
  OptFunction f(
      [&, this](const std::vector<double> &x, std::vector<double> &dx) {

        auto evaled = kernel->operator()(x);
        auto kernels = uniqueTerms->observe(evaled);
        auto nEnergyKernels = kernels.size();

        // now get circuits for gradients of unique terms
        int nGradKernels = 0;
        if (optimizer->isGradientBased()) {
          auto parameterShift = xacc::getGradient("parameter-shift");
          parameterShift->initialize(
              {{"observable", uniqueTerms}, {"shift-scalar", 0.5}});
          auto gradKernels = parameterShift->getGradientExecutions(kernel, x);

          nGradKernels = gradKernels.size();
          kernels.insert(kernels.end(), gradKernels.begin(), gradKernels.end());
        }

        // excecute all circuits
        auto tmpBuffer = qalloc(buffer->size());
        accelerator->execute(tmpBuffer, kernels);

        nTotalCircuits += kernels.size();
        for (auto &op : uniqueTerms->getNonIdentitySubTerms()) {
          nTotalMeasurements += std::dynamic_pointer_cast<PauliOperator>(op)
                                    ->getTerms()
                                    .begin()
                                    ->second.ops()
                                    .size();
        }

        // keep track of buffers
        auto buffers = tmpBuffer->getChildren();
        auto energyBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(
            buffers.begin(), buffers.begin() + nEnergyKernels);
        auto gradBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(
            buffers.begin() + nEnergyKernels,
            buffers.begin() + nEnergyKernels + nGradKernels);

        // compute moments and spectrum
        auto moments = computeMoments(momentOperators, energyBuffers);
        auto spectrum = PDS(moments);

        // get ground state energy
        auto energy = *std::min_element(spectrum.begin(), spectrum.end());
        if (energy < optEnergy) {
          optEnergy = energy;
          optSpectrum = spectrum;
        }

        // get denominator and derivative of the polynomial in eq. 17
        Eigen::VectorXd polynomialDerivative =
            getPolynomialDerivative(X, energy);

        // get moments gradients
        auto nParams = x.size();
        auto nTerms = uniqueTerms->getNonIdentitySubTerms().size();
        auto momentsGradients = computeMomentsGradients(
            momentOperators, gradBuffers, nTerms, nParams);

        // get the ground state energy gradient
        Eigen::VectorXd dE = Eigen::VectorXd::Zero(nParams);
        for (int i = 0; i < nParams; i++) {

          // get gradient wrt parameter
          auto gradMomentsParam = momentsGradients[i];

          // get dM
          Eigen::MatrixXd dM(order, order);
          dM = computeMatrixM(gradMomentsParam, true);
          // M(order, order) = <trial|trial> = 1
          // so it has no angle dependence, thus its derivative is 0

          // get dY
          Eigen::VectorXd dY(order);
          dY = computeVectorY(gradMomentsParam);

          // get dX
          Eigen::VectorXd dX(order);
          dX = M.colPivHouseholderQr().solve(-dY - dM * X).transpose();

          // get dE
          dE(i) = polynomialDerivative.dot(dX);
        }

        dx = std::vector<double>(dE.data(), dE.data() + dE.size());
        ss << "PDS-VQS(" << order << "): " << energy << "\t dx: " << dx << "\n";
        xacc::info(ss.str());
        ss.str(std::string());

        return energy;
      },
      kernel->nVariables());

  if (!adapt) {

    auto result = optimizer->optimize(f);
    buffer->addExtraInfo("spectrum", optSpectrum);
    buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
    buffer->addExtraInfo("opt-params", ExtraInfo(result.second));

    ss << "Total number of circuit implementations: " << nTotalCircuits;
    xacc::info(ss.str());
    ss.str(std::string());

    ss << "Total number of measurements: " << nTotalMeasurements;
    xacc::info(ss.str());
    ss.str(std::string());

  } else {

    // get operator pool
    auto pool = xacc::getService<OperatorPool>(poolName);
    pool->optionalParameters({{"n-electrons", nElectrons}});
    auto operators = pool->generate(buffer->size());

    // get operator for moments gradients
    std::vector<std::shared_ptr<Observable>> commutators;
    for (auto &op : operators) {
      for (auto &momentOperator : momentOperators) {
        commutators.push_back(momentOperator->commutator(op));
      }
    }

    // get unique terms for gradients
    std::vector<std::shared_ptr<Observable>> uniqueGradientStrings(commutators);
    uniqueGradientStrings.insert(uniqueGradientStrings.end(),
                                 momentOperators.begin(),
                                 momentOperators.end());
    auto uniqueStringsPtr = getUniqueTerms(uniqueGradientStrings);
    auto nUniqueStrings = uniqueStringsPtr->getNonIdentitySubTerms().size();

    auto nMeasurements = 0;
    for (auto &op : uniqueStringsPtr->getNonIdentitySubTerms()) {
      nMeasurements += std::dynamic_pointer_cast<PauliOperator>(op)
                           ->getTerms()
                           .begin()
                           ->second.ops()
                           .size();
    }

    // start adapt
    std::vector<double> x;
    double oldEnergy = 0.0;
    for (int iter = 0; iter < maxIterations; iter++) {

      // observe unique strings for current ansatz
      auto evaled = kernel->operator()(x);
      auto observedKernels = uniqueStringsPtr->observe(evaled);
      nTotalMeasurements += nMeasurements;
      nTotalCircuits += observedKernels.size();

      // execute all unique terms
      auto tmpBuffer = qalloc(buffer->size());
      accelerator->execute(tmpBuffer, observedKernels);
      auto buffers = tmpBuffer->getChildren();

      // solve PDS for the current ansatz
      auto currentMoments = computeMoments(momentOperators, buffers);
      auto spectrum = PDS(currentMoments);
      auto energy = *std::min_element(spectrum.begin(), spectrum.end());
      Eigen::VectorXd polynomialDerivative = getPolynomialDerivative(X, energy);

      // get gradients
      auto begin = commutators.begin();
      int maxGradientIdx = 0;
      double maxGradient = 0.0, gradientNorm = 0.0;
      for (int i = 0; i < operators.size(); i++) {

        auto momentsGradientOperator = std::vector<std::shared_ptr<Observable>>(
            begin, begin + 2 * order - 1);
        auto momentsGradient = computeMoments(momentsGradientOperator, buffers);
        begin += 2 * order - 1;

        // get matrices/vectors for commutators
        Eigen::MatrixXd gradientM = computeMatrixM(momentsGradient, true);
        Eigen::VectorXd gradientY = computeVectorY(momentsGradient);
        Eigen::VectorXd gradientX = M.colPivHouseholderQr()
                                        .solve(-gradientY - gradientM * X)
                                        .transpose();
        auto gradient = polynomialDerivative.dot(gradientX);

        // print gradient elements above threshold
        if (abs(gradient) > printThreshold) {
          ss << "dE/dt" << i << " = " << gradient;
          xacc::info(ss.str());
          ss.str(std::string());
        }

        // update max gradient
        if (abs(gradient) > abs(maxGradient)) {
          maxGradientIdx = i;
          maxGradient = gradient;
        }

        gradientNorm += gradient * gradient;
      }

      ss << "Max gradient component: dE/dt " << maxGradientIdx
         << " = " << maxGradient << " a.u.";
      xacc::info(ss.str());
      ss.str(std::string());

      gradientNorm = std::sqrt(gradientNorm);
      ss << "Norm of gradient vector: " << gradientNorm << " a.u.";
      xacc::info(ss.str());
      ss.str(std::string());

      if (gradientNorm < _adaptThreshold) {
        ss << "Converged!";
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Total number of circuit implementations: " << nTotalCircuits;
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Total number of measurements: " << nTotalMeasurements;
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Depth of final circuit: " << kernel->depth();
        xacc::info(ss.str());
        ss.str(std::string());

        return;

      } else {

        auto maxGradientGate =
            pool->getOperatorInstructions(maxGradientIdx, kernel->nVariables());
        kernel->addVariable("x" + std::to_string(kernel->nVariables()));
        for (auto &inst : maxGradientGate->getInstructions()) {
          kernel->addInstruction(inst);
        }
        x.push_back(0.0);

        optimizer->appendOption("initial-parameters", x);
        auto result = optimizer->optimize(f);

        auto newEnergy = result.first;
        oldEnergy = newEnergy;
        ss << "Energy at ADAPT iteration " << iter + 1 << ": " << newEnergy;
        xacc::info(ss.str());
        ss.str(std::string());

        x = result.second;
        ss << "Parameters at ADAPT iteration " << iter + 1 << ": ";
        for (auto param : x) {
          ss << param << " ";
        }
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Circuit depth at ADAPT iteration " << iter + 1 << ": "
           << kernel->depth();
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Number of circuit implementations at ADAPT iteration "
           << iter + 1 << ": " << nTotalCircuits;
        xacc::info(ss.str());
        ss.str(std::string());

        ss << "Number of measurements at ADAPT iteration " << iter + 1 << ": "
           << nTotalMeasurements;
        xacc::info(ss.str());
        ss.str(std::string());

        buffer->addExtraInfo("spectrum", optSpectrum);
        buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
        buffer->addExtraInfo("opt-params", ExtraInfo(result.second));
      }
    }
  }

  return;
}

// Compute energy from CMX
// J. Phys. A: Math. Gen. 17, 625 (1984)
// and Int. J. Mod. Phys. B 9, 2899 (1995)
std::vector<double> PDS_VQS::PDS(const std::vector<double> &moments) const {

  M = computeMatrixM(moments, false);
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

Eigen::MatrixXd PDS_VQS::computeMatrixM(const std::vector<double> &moments,
                                        const bool gradient) const {

  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(order, order);
  for (int i = 0; i < order; i++) {
    for (int j = i; j < order; j++) {
      if (i + j == 2 * (order - 1)) {
        ret(i, j) = gradient ? 0.0 : 1.0;
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

Eigen::VectorXd PDS_VQS::getPolynomialDerivative(const Eigen::VectorXd X,
                                                 const double energy) const {

  // get denominator and derivative of the polynomial in eq. 17
  auto denominator = -order * std::pow(energy, order - 1);
  Eigen::VectorXd polynomialDerivative = Eigen::VectorXd::Zero(order);
  polynomialDerivative(order - 1) = 1.0;
  for (auto i = 0; i < order - 1; i++) {
    denominator -= (order - i - 1) * X(i) * std::pow(energy, order - i - 2);
    polynomialDerivative(i) = std::pow(energy, order - i - 1);
  }
  polynomialDerivative /= denominator;

  return polynomialDerivative;
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
           buffers[termIdx + (2 * paramIdx + 1) * nTerms]
               ->getExpectationValueZ()) /
          2.0;
      paramUniqueTermGrad[termIdx] = {termName, termGrad};
    }

    std::vector<double> momentsGradientsParams(momentOperators.size());
    for (int m = 0; m < 2 * order - 1; m++) {

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

} // namespace algorithm
} // namespace xacc
