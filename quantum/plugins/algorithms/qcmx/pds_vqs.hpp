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
// This implements the Quantum Connected Moment Expansion (PDS_VQS) algorithm
// This implementation follows arXiv:2009.05709v2
//
#ifndef XACC_ALGORITHM_PDS_VQS_HPP_
#define XACC_ALGORITHM_PDS_VQS_HPP_

#include "Algorithm.hpp"
#include <Eigen/Dense>
#include "Observable.hpp"

namespace xacc {
namespace algorithm {

class PDS_VQS : public Algorithm {

protected:
  std::shared_ptr<Observable> observable;
  Accelerator *accelerator;
  CompositeInstruction *kernel;

  // CMX order, also K in the paper
  int maxOrder;
  // threshold below which we ignore measurement
  double threshold = 0.0;
  // x is the vector of parameters if the provided ansatz is parameterized
  std::vector<double> x;
  // spectrum is only for PDS
  mutable std::vector<double> spectrum;
  // store the matrices for gradients
  mutable Eigen::MatrixXd M;
  mutable Eigen::MatrixXd X;
  mutable Eigen::VectorXd Y;
  //mutable Eigen::VectorXd polynomial;

  double measureOperator(const std::shared_ptr<Observable> obs,
                         const int bufferSize) const;

  // Compute energy from PDS CMX
  double PDS(const std::vector<double> &moments, const int order) const;

  // get unique Pauli words
  std::shared_ptr<Observable>
  getUniqueTerms(const std::vector<std::shared_ptr<Observable>>) const;

  // get moments from measured unique Pauli words
  std::vector<double>
  computeMoments(const std::vector<std::shared_ptr<Observable>>,
             const std::vector<std::shared_ptr<AcceleratorBuffer>>) const;

  // compute matrix M
  Eigen::MatrixXd computeMatrixM(const std::vector<double> &) const;
  // compute vector Y
  Eigen::VectorXd computeVectorY(const std::vector<double> &) const;

  Eigen::MatrixXd computeMomentsGradients(
    const std::shared_ptr<Observable>,
    const std::vector<std::shared_ptr<Observable>>,
    const std::vector<std::shared_ptr<AcceleratorBuffer>> , const int ) const ;

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