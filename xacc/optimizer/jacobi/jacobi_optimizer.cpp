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
//
// This implements several optimization strategies
// presented in in arXiv:1904.03206
//
// This paper represents Ry gates as
// Ry = exp(-i * x * Y)
// while in XACC
// Ry = exp((-i * x * Y)/2)
// So all equations used in this implementation were rewritten
// in order to comply with our convention

#include "jacobi_optimizer.hpp"
#include "CompositeInstruction.hpp"
#include "xacc.hpp"
#include <unsupported/Eigen/KroneckerProduct>
#include <cmath>
#include <functional>
#include <numeric>

namespace xacc {

const std::string JacobiOptimizer::get_algorithm() const {

  std::string jacobi_opt_name = "jacobi-1";
  if (options.stringExists("jacobi-optimizer")) {
    jacobi_opt_name = options.getString("jacobi-optimizer");
  }

  return jacobi_opt_name;
}

OptResult JacobiOptimizer::optimize(OptFunction &function) {

  if (!options.pointerLikeExists<CompositeInstruction>("ansatz")) {
    xacc::error("Jacobi optimizer requires the ansatz");
  }
  ansatz = options.getPointerLike<CompositeInstruction>("ansatz");
  getVariableCoefficients();

  // default strategy
  std::string jacobi_opt_name = "jacobi-1";
  if (options.stringExists("jacobi-optimizer")) {
    jacobi_opt_name = options.getString("jacobi-optimizer");
  }

  // number of parameterized gates
  nParams = function.dimensions();
  x.resize(nParams, 0.0);
  grad.resize(nParams, 0.0);
  gateCoefficients.resize(nParams, 0.0);

  // check for acceleration
  if (options.stringExists("acceleration")) {
    acceleration = options.getString("acceleration");

    // set up DIIS
    if (options.keyExists<int>("n-vector-history")) {
      maxMetricHistory = options.get<int>("n-vector-history");
    }
    metricHistory.resize(nParams, 1);
  }

  // T^-1 Eq.10
  transferMatrixInverse1 << 1.0, 0.0, 1.0, -1.0, 2.0, -1.0, -1.0, 0.0, 1.0;
  transferMatrixInverse1 /= 2.0;

  if (jacobi_opt_name == "jacobi-1") {

    if (options.keyExists<std::vector<int>>("cluster-sequence")) {
      clusterSeq1 = options.get<std::vector<int>>("cluster-sequence");
    } else {
      clusterSeq1.resize(nParams);
      std::iota(clusterSeq1.begin(), clusterSeq1.end(), 0);
    }

    return jacobi1(function);
  }

  if (jacobi_opt_name == "jacobi-2") {

    twoGateGridExpValues.resize(9);

    // T^-1 for Jacobi-2 
    transferMatrixInverse2 =
        Eigen::kroneckerProduct(transferMatrixInverse1, transferMatrixInverse1);

    if (options.keyExists<std::vector<std::pair<int, int>>>(
            "cluster-sequence")) {
      clusterSeq2 =
          options.get<std::vector<std::pair<int, int>>>("cluster-sequence");

    } else {

      for (int i = 1; i < nParams; i++) {
        for (int j = 0; j < i; j++) {
          clusterSeq2.push_back(std::make_pair(j, i));
        }
      }

    }

    return jacobi2(function);
  }
}

// This implements the Jacobi-1 strategy
OptResult JacobiOptimizer::jacobi1(OptFunction &function) {

  std::vector<double> newX(nParams), oldMetric(nParams), newMetric(nParams);
  double newEnergy, oldEnergy;

  // Compute <H(0.0,..., 0.0)> only once
  oneGateGridExpValues(1) = function(x, grad);

  _maxiter = 2;
  // loop until norm of error/gradient vector is small
  for (int iter = 0; iter < _maxiter; iter++) {

    // error vector norm
    double norm = 0.0;

    // compute the objective function for the grid points
    // for each parameter in the clusterSeq
    for (int i : clusterSeq1) {

      // seqX stores the new x for each element in x
      auto seqX = x;

      // Get optimal value for parameter i
      // from single gate tomography
      seqX[i] = oneGateTomography(i, x, function, newMetric);

      // Compute energy with new parameter and store it
      newEnergy = function(seqX, grad);
      newX[i] = seqX[i];
    }

    // compute norm of the error metric
    for (int i = 0; i < nParams; i++) {
      norm += std::pow(newMetric[i] - oldMetric[i], 2);
    }
    norm = std::sqrt(norm);

    // if norm is small, stop
    if (norm < _threshold)
      break;

    // Use acceleration (Anderson or Pulay)
    // add error/gradient vector to history
    // delete oldest vector if necessary
    // pass history to the DIIS function
    if (iter > 1 && !acceleration.empty()) {

      Eigen::Map<Eigen::VectorXd> eigenMetric(newMetric.data(),
                                              newMetric.size());

      if (metricHistory.cols() < maxMetricHistory) {

        metricHistory.conservativeResize(Eigen::NoChange,
                                         metricHistory.cols() + 1);
        metricHistory.col(metricHistory.cols() - 1) = eigenMetric;

      } else {

        metricHistory.block(0, 0, metricHistory.rows() - 1,
                            metricHistory.cols() - 2) =
            metricHistory.block(0, 1, metricHistory.rows() - 1,
                                metricHistory.cols() - 1);
        metricHistory.col(metricHistory.cols() - 1) = eigenMetric;
      }

    } else {
      x = newX;
    }
  }

  OptResult result(newEnergy, newX);

  return result;
}

// this is not an analytical solution
// this is one way by which we can estimate the parameter update
// for simultaneous optimization of parameters for two gates
OptResult JacobiOptimizer::jacobi2(OptFunction &function) {

  std::vector<double> newX(nParams), oldMetric(nParams),
      newMetric(nParams), energies(nParams, 0.0);
  double zeroExpValue, newEnergy;

  // Compute <H(0.0,..., 0.0)> only once
  zeroExpValue = function(x, grad);
  twoGateGridExpValues(4) = zeroExpValue;
  oneGateGridExpValues(1) = zeroExpValue;

  _maxiter = 1;

  for (int iter = 0; iter < _maxiter; iter++) {

    // error vector norm
    double norm = 0.0;

    // compute the objective function for the grid points
    // for each parameter in the clusterSeq
    for (auto pair : clusterSeq2) {

      auto seqX = x;
      double oldEnergy = 0.0;
      newEnergy = zeroExpValue;

      // update B
      seqX[pair.second] =
          oneGateTomography(pair.second, x, function, newMetric);

      
      while (fabs(newEnergy - oldEnergy) > energyThreshold) {
      //for (int i =0;i<1;i++){

        oldEnergy = newEnergy;
        // update A with new B
        seqX[pair.first] = twoGateTomography(pair.first, pair.second, seqX,
                                             function, newMetric);
        // update B with new A
        seqX[pair.second] = twoGateTomography(pair.second, pair.first, seqX,
                                              function, newMetric);

        newEnergy = function(seqX, grad);
      }

      if (newEnergy < energies[pair.first]){
        newX[pair.first] = seqX[pair.first];
        energies[pair.first] = newEnergy;
      }

      if (newEnergy < energies[pair.second]){
        newX[pair.second] = seqX[pair.second];
        energies[pair.second] = newEnergy;
      }
    }
  }

  OptResult result(newEnergy, newX);
  return result;
}

std::vector<double> JacobiOptimizer::DIIS(Eigen::MatrixXd metric) {

  // Construct overlap of error vectors
  auto nCols = metric.cols();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nCols + 1, nCols + 1);
  B.block(0, 0, nCols - 1, nCols - 1) = metric.transpose() * metric;
  for (int i = 0; i < nCols; i++) {
    B(i, nCols + 1) = -1.0;
    B(nCols + 1, i) = -1.0;
  }

  // RHS
  Eigen::VectorXd x = Eigen::VectorXd::Zero(nCols + 1);
  x(nParams + 1) = -1.0;

  // Solve Bc = x
  Eigen::VectorXd extrapCoeffs = B.colPivHouseholderQr().solve(x);

  // parameter extrapolation
  Eigen::VectorXd parameterUpdate = Eigen::VectorXd::Zero(nCols);
  for (int i = 0; i < nCols; i++) {
    parameterUpdate += extrapCoeffs(i) * metric.col(i);
  }

  std::vector<double> newX(parameterUpdate.data(),
                           parameterUpdate.data() + parameterUpdate.size());

  return newX;
}

// returns optimal value of gate parameter from tomography
double JacobiOptimizer::oneGateTomography(const int paramIdx,
                                          const std::vector<double> x,
                                          OptFunction &function,
                                          std::vector<double> &metric) {

  // loop over grid points and compute the objective function Eq.8
  for (int point : {0, 2}) {

    auto xForGrid = x;
    xForGrid[paramIdx] = gridPoints[point] * gateCoefficients[paramIdx];
    oneGateGridExpValues(point) = function(xForGrid, grad);
  }

  // Compute the tomography parameters Eq.9
  Eigen::Vector3d tomoParams = transferMatrixInverse1 * oneGateGridExpValues;

  double optimalParam =
      atan2(-tomoParams(2), -tomoParams(1));// * gateCoefficients[paramIdx];

  if (acceleration == "pulay") {
    metric[paramIdx] =
        -tomoParams(1) * sin(x[paramIdx]) + tomoParams(2) * cos(x[paramIdx]);
  } else {
    metric[paramIdx] = optimalParam;
  }
  return optimalParam;
}

// returns optimal value of gate parameter from tomography
// this will return the update for parameter paramIdx1 for a given
// value of paramIdx2
double JacobiOptimizer::twoGateTomography(const int paramIdx1,
                                          const int paramIdx2,
                                          const std::vector<double> x,
                                          OptFunction &function,
                                          std::vector<double> &metric) {

  // loop over grid points and compute the objective function Eq.8
  int point = 0;
  for (int point1 : {0, 1, 2}) {

    auto xForGrid = x;
    xForGrid[paramIdx1] = gridPoints[point1] * gateCoefficients[paramIdx1];
    for (int point2 : {0, 1, 2}) {
      if ((point1 == 1) && (point2 == 1)) {
        point++;
        continue;
      }

      xForGrid[paramIdx2] = gridPoints[point2] * gateCoefficients[paramIdx2];
      twoGateGridExpValues(point) = function(xForGrid, grad);

      point++;
    }
  }

  // Compute the tomography parameters Eq.9
  Eigen::VectorXd tomoParams = transferMatrixInverse2 * twoGateGridExpValues;

  double numerator = tomoParams(6) + tomoParams(7) * cos(x[paramIdx2]) +
                     tomoParams(8) * sin(x[paramIdx2]);
  double denominator = tomoParams(3) + tomoParams(4) * cos(x[paramIdx2]) +
                       tomoParams(5) * sin(x[paramIdx2]);

  double optimalParam =
      atan2(-numerator, -denominator) * gateCoefficients[paramIdx1];

  return optimalParam;
}

  // This returns the coefficients of the Pauli strings
  // in parameterized gates
  void JacobiOptimizer::getVariableCoefficients() {

    // get all the variables and sort them
    auto vars = ansatz->getVariables();
    auto kernelInst = ansatz->getInstructions();

    // loop over variables
    // look for them in the instructions
    // once found, retrieve the coefficient
    for (auto &var : vars) {

      // this bool works as a continue for the outer loop
      bool hasFound = false;
      for (int i = 0; !hasFound && i < kernelInst.size(); i++) {

        if (kernelInst[i]->isParameterized() &&
            kernelInst[i]->getParameter(0).isVariable()) {

          auto paramStr = kernelInst[i]->getParameter(0).as<std::string>();
          if (var == paramStr) {

            gateCoefficients.push_back(1.0);
            hasFound = true;

          } else if (paramStr.find(var) != std::string::npos) {

            auto coeff = std::stod(paramStr.substr(0, paramStr.find("*")));
            gateCoefficients.push_back(fabs(coeff));
            hasFound = true;
          }
        }
      }
    }

  }

} // namespace xacc
