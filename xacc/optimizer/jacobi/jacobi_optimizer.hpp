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
 *******************************************************************************/
#ifndef XACC_JACOBI_OPTIMIZER_HPP_
#define XACC_JACOBI_OPTIMIZER_HPP_

#include <memory>
#include <type_traits>
#include <utility>
#include <Eigen/Dense>
#include "CompositeInstruction.hpp"
#include "Optimizer.hpp"
#include "xacc.hpp"

namespace xacc {

class JacobiOptimizer : public Optimizer {
public:
  OptResult optimize(OptFunction &function) override;
  const bool isGradientBased() const override {return true;};
  virtual const std::string get_algorithm() const;

  const std::string name() const override { return "jacobi"; }
  const std::string description() const override { return ""; }

protected:

  std::string acceleration = "";
  double _threshold = 1.0e-3, energyThreshold = 1.0e-6;
  int nParams, _maxiter = 50, maxMetricHistory = 5;
  std::vector<int> clusterSeq1;
  std::vector<std::pair<int,int>> clusterSeq2;
  std::vector<double> gridPoints{-xacc::constants::pi / 2.0, 0.0, xacc::constants::pi / 2.0};
  std::vector<double> x, grad, errorVector, gateCoefficients;
  Eigen::MatrixXd metricHistory, transferMatrixInverse2;
  Eigen::Matrix3d transferMatrixInverse1;
  Eigen::Vector3d oneGateGridExpValues;
  Eigen::VectorXd twoGateGridExpValues;
  CompositeInstruction* ansatz;

  OptResult jacobi1(OptFunction &function);
  
  OptResult jacobi2(OptFunction &function);

  std::vector<double> DIIS(Eigen::MatrixXd metric);

  double oneGateTomography(const int paramIdx, const std::vector<double> x, OptFunction &function, std::vector<double>& metric);

  double twoGateTomography(const int paramIdx1, const int paramIdx2, const std::vector<double> x, OptFunction &function, std::vector<double>& metric);

  void getVariableCoefficients();

};
} // namespace xacc
#endif
