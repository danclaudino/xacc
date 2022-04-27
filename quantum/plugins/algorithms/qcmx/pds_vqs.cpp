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
#include "Utils.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "PauliOperator.hpp"

#include <algorithm>
#include <complex>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

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

  accelerator = parameters.getPointerLike<Accelerator>("accelerator");
  order = parameters.get<int>("cmx-order");
  kernel = parameters.getPointerLike<CompositeInstruction>("ansatz");

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

  // in case the ansatz is parameterized
  if (parameters.keyExists<std::vector<double>>("parameters")) {
    x = parameters.get<std::vector<double>>("parameters");
    std::reverse(x.begin(), x.end());
  }

  return true;
}

const std::vector<std::string> PDS_VQS::requiredParameters() const {
  return {"observable", "accelerator", "cmx-order", "ansatz"};
}

void PDS_VQS::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  // First gather the operators for all required moments
  std::vector<std::shared_ptr<Observable>> momentOperators(2 * order - 1);
  auto momentOperator = PauliOperator("I");
  for (int i = 0; i < 2 * order - 1; i++) {
    momentOperator *= (*std::dynamic_pointer_cast<PauliOperator>(observable));
    momentOperators.push_back(std::make_shared<PauliOperator>(momentOperator));
  }

  // now get the unique terms in all moments operators
  std::vector<double> x{7 * xacc::constants::pi / 16.0, xacc::constants::pi, 0, 0};
  auto uniqueTerms = getUniqueTerms(momentOperators);
  auto evaled = kernel->operator()(x);
  auto kernels = uniqueTerms->observe(evaled);
  auto nEnergyKernels = kernels.size();
  // now get executions for gradients of unique terms
  auto parameterShift = xacc::getGradient("parameter-shift");
  parameterShift->initialize({{"observable", uniqueTerms}});
  auto gradKernels = parameterShift->getGradientExecutions(xacc::as_shared_ptr(kernel), x);
  kernels.insert(kernels.end(), gradKernels.begin(), gradKernels.end());

  // excecute
  accelerator->execute(buffer, kernels);

  // keep track of buffers
  auto buffers = buffer->getChildren();
  auto energyBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(buffers.begin(), buffers.begin() + nEnergyKernels);
  auto gradBuffers = std::vector<std::shared_ptr<AcceleratorBuffer>>(buffers.begin() + nEnergyKernels, buffers.end());

  // compute moments and spectrum
  auto moments = computeMoments(momentOperators, energyBuffers);
  auto spectrum = PDS(moments);

  // get ground state energy
  auto energy = *std::min_element(spectrum.begin(), spectrum.end());

  std::cout << energy << "\n";
  
  exit(0);

  // get vector that dot with dX/dx in eq. 17
  auto denominator = order * std::pow(energy, order - 1);
  Eigen::VectorXd polynomial(order - 1);
  polynomial(order - 1) = 1.0;
  for (auto i = 1; i < order; i++) {
    denominator += (order - i) * X(i) * std::pow(energy, order - i - 1);
    polynomial(i - 1) = std::pow(energy, order - i - 1);
  }
  polynomial /= denominator;

  // get moments gradients
  // it is a matrix with params x moments
  Eigen::MatrixXd momentsGradients(x.size(), momentOperators.size());
  momentsGradients = computeMomentsGradients(uniqueTerms, momentOperators, gradBuffers, x.size());

  // get the ground state energy gradient
  Eigen::VectorXd dE(x.size());
  for (int i = 0; i < x.size(); i++) {

    // get dM
    std::vector<double> gradMomentsParam(momentsGradients.row(i).data(), momentsGradients.row(i).data() + momentsGradients.row(i).size());
    Eigen::MatrixXd dM(order, order);
    dM = computeMatrixM(gradMomentsParam);

    // get dY
    Eigen::VectorXd dY(order);
    dY = computeVectorY(gradMomentsParam);

    // get dX
    Eigen::VectorXd dX(order);
    dX = M.colPivHouseholderQr().solve(-dY + dM * X).transpose();

    // get dE
    dE(i) = polynomial.dot(dX);

  }

  // apply metric




  // add energy to the buffer
  buffer->addExtraInfo("spectrum", spectrum);
  return;
}

// Compute energy from CMX
// J. Phys. A: Math. Gen. 17, 625 (1984)
// and Int. J. Mod. Phys. B 9, 2899 (1995)
std::vector<double> PDS_VQS::PDS(const std::vector<double> &moments) const {

  M = computeMatrixM(moments);
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
  for (auto &momentOp : momentOps) {
    for (auto &subTerm : momentOp->getNonIdentitySubTerms()) {

      if (!xacc::container::contains(uniqueTermsStrings, subTerm->toString())) {

        if (std::fabs(subTerm->coefficient().real()) > threshold) {
          uniqueTermsStrings.push_back(subTerm->toString());
        }
      }
    }
  }

  auto uniqueTermsPtr = std::make_shared<PauliOperator>();
  for (auto &pauliStr : uniqueTermsStrings) {
    uniqueTermsPtr->operator+=(*std::make_shared<PauliOperator>(pauliStr));
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
          expval += std::real(subTerm->coefficient() * buffer->getExpectationValueZ());
          break;
        }
      }
    }
    moments.push_back(expval);
  }
  return moments;
}

Eigen::MatrixXd PDS_VQS::computeMatrixM(const std::vector<double> &moments) const {

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

Eigen::VectorXd PDS_VQS::computeVectorY(const std::vector<double> &moments) const {

  Eigen::VectorXd ret = Eigen::VectorXd::Zero(order);
  for (int i = 0; i < order; i++) {
    ret(i) = moments[2 * order - i - 2];
  }

  return ret;
}

Eigen::MatrixXd PDS_VQS::computeMomentsGradients(
    const std::shared_ptr<Observable> uniqueTerms,
    const std::vector<std::shared_ptr<Observable>> momentOperators,
    const std::vector<std::shared_ptr<AcceleratorBuffer>> buffers, const int nParams) const {

    // loop over all parameters
    auto terms = uniqueTerms->getNonIdentitySubTerms();
    auto nTerms = terms.size();
    //std::vector<std::vector<std::pair<std::string, double>>> uniqueTermsGradPerParam(nParams);
    Eigen::MatrixXd momentsGradients = Eigen::MatrixXd::Zero(nParams, momentOperators.size());
    //std::vector<std::vector<double>> momentsGradients(nParams);
    for (int paramIdx = 0; paramIdx < nParams; paramIdx++) {

      std::vector<std::pair<std::string, double>> paramUniqueTermGrad(nTerms);
      for (int termIdx = 0; termIdx < nTerms; termIdx++) {

        auto termName = buffers[termIdx + 2 * nParams]->name();
          
        auto termGrad = (buffers[termIdx + 2 * nParams]->getExpectationValueZ() - buffers[termIdx + 2 * nParams + nTerms]->getExpectationValueZ()) / 2.0;

        paramUniqueTermGrad[termIdx] = {termName, termGrad};
      }

      //uniqueTermsGradPerParam[paramIdx] = paramUniqueTermGrad;
      std::vector<double> momentsGradientsParams(momentOperators.size());
      int momentCounter = 0;
      for (auto & momentOp : momentOperators) {

        double momentGrad = 0.0;
        for (auto & subTerm : momentOp->getNonIdentitySubTerms()) {

          auto termName = std::dynamic_pointer_cast<PauliOperator>(subTerm)->begin()->first;
          auto coeff = subTerm->coefficient().real();

          for (auto & pair : paramUniqueTermGrad) {

            if (termName == pair.first) {
              momentGrad += coeff * pair.second;
              break;
            }

          }

        }
       // momentsGradientsParams.push_back(momentGrad);
       momentsGradients(paramIdx, momentCounter++) = momentGrad;

      }
      //momentsGradients[paramIdx] = momentsGradientsParams;
    }



/*
  std::vector<double> moments;
  for (int i = 0; i < 2 * maxOrder - 1; i++) {
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
          expval += std::real(subTerm->coefficient() * buffer->getExpectationValueZ());
          break;
        }
      }
    }
    moments.push_back(expval);
  }
  */
  return momentsGradients;


    }

} // namespace algorithm
} // namespace xacc
