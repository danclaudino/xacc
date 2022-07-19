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
// This implements the PDS_VQS algorithm
// This implementation follows https://doi.org/10.22331/q-2021-06-10-473
//
#ifndef XACC_ALGORITHM_PDS_VQS_HPP_
#define XACC_ALGORITHM_PDS_VQS_HPP_

#include "Algorithm.hpp"
#include <Eigen/Dense>

namespace xacc {
namespace algorithm {

class PDS_VQS : public Algorithm {

protected:
  std::shared_ptr<Observable> observable;
  Accelerator *accelerator;
  CompositeInstruction *initialState;
  Optimizer* optimizer;

  bool adapt = false;
  std::string poolName = "";
  int nElectrons, maxIterations = 20;
  double _adaptThreshold = 1.0e-2;

  // CMX order, also K in the paper
  int order;
  // threshold below which we ignore measurement
  double printThreshold = 1.0e-5;

  // store the matrices for gradients
  mutable Eigen::MatrixXd M;
  mutable Eigen::VectorXd X;

  // Compute energy from PDS CMX
  std::vector<double> PDS(const std::vector<double> &moments) const;

  // get unique Pauli words
  std::shared_ptr<Observable>
  getUniqueTerms(const std::vector<std::shared_ptr<Observable>>) const;

  // get moments from measured unique Pauli words
  std::vector<double>
  computeMoments(const std::vector<std::shared_ptr<Observable>>,
                 const std::vector<std::shared_ptr<AcceleratorBuffer>>) const;

  // compute matrix M
  Eigen::MatrixXd computeMatrixM(const std::vector<double> &, const bool) const;
  // compute vector Y
  Eigen::VectorXd computeVectorY(const std::vector<double> &) const;

  // get derivative of the polynomial
  Eigen::VectorXd getPolynomialDerivative(const Eigen::VectorXd, const double) const;

  //Eigen::MatrixXd
  std::vector<std::vector<double>>
  computeMomentsGradients(const std::vector<std::shared_ptr<Observable>>,
                          const std::vector<std::shared_ptr<AcceleratorBuffer>>,
                          const int, const int) const;

  Eigen::MatrixXd getMetricMatrix(const std::vector<double> &, const int) const;

public:
  bool initialize(const HeterogeneousMap &parameters) override;
  const std::vector<std::string> requiredParameters() const override;
  void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
  const std::string name() const override { return "pds-vqs"; }
  const std::string description() const override { return ""; }

  DEFINE_ALGORITHM_CLONE(PDS_VQS)
};

//
} // namespace algorithm
} // namespace xacc

#endif