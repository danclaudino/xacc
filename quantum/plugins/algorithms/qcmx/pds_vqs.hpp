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
#include "Optimizer.hpp"
#include "PauliOperator.hpp"
#include <Eigen/Dense>
#include <vector>
#include "Eigen/src/Core/Map.h"
#include "Observable.hpp"

namespace xacc {
namespace algorithm {

// this is from Thien's QNG implementation
struct ParametrizedCircuitLayer {
  // Gates that precede the layer
  std::vector<InstPtr> preOps;
  // Parametrized operators in the layer
  std::vector<InstPtr> ops;
  // Corresponding optimization parameter indices
  std::vector<size_t> paramInds;
  // Gates that succeed the layer
  std::vector<InstPtr> postOps;
  // Tensor matrix terms:
  std::vector<xacc::quantum::PauliOperator> kiTerms;
  std::vector<xacc::quantum::PauliOperator> kikjTerms;
  // Partition the circuit into layers.
  static std::vector<ParametrizedCircuitLayer> toParametrizedLayers(
      const std::shared_ptr<xacc::CompositeInstruction> &in_circuit);
};
using ObservedKernels =
    std::vector<std::shared_ptr<xacc::CompositeInstruction>>;

class PDS_VQS : public Algorithm {

protected:
  std::shared_ptr<Observable> observable;
  Accelerator *accelerator;
  CompositeInstruction *kernel;
  Optimizer* optimizer;

  // CMX order, also K in the paper
  int order, nRoots = 1;
  // threshold below which we ignore measurement
  double threshold = 0.0, step = 0.05;

  // store the matrices for gradients
  mutable Eigen::MatrixXd M;
  mutable Eigen::MatrixXd X;
  //mutable Eigen::VectorXd Y;

  std::string metric = "GD";
  std::vector<std::string> IMPLD_METRICS{"GD", "NGD", "ITE"};

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
  Eigen::MatrixXd computeMatrixM(const std::vector<double> &) const;
  // compute vector Y
  Eigen::VectorXd computeVectorY(const std::vector<double> &) const;

  //Eigen::MatrixXd
  std::vector<std::vector<double>>
  computeMomentsGradients(const std::vector<std::shared_ptr<Observable>>,
                          const std::vector<std::shared_ptr<AcceleratorBuffer>>,
                          const int, const int) const;

  Eigen::MatrixXd getMetricMatrix(const std::vector<double> &, const int) const;

  ObservedKernels
  getMetricCircuits(std::shared_ptr<xacc::CompositeInstruction> in_circuit,
                    const std::vector<double> &in_x) const;

  ObservedKernels
  constructMetricTensorSubCircuit(ParametrizedCircuitLayer &io_layer,
                                  const std::vector<std::string> &in_varNames,
                                  const std::vector<double> &in_varVals) const;

  Eigen::MatrixXd constructMetricTensorMatrix(
      const std::vector<std::shared_ptr<xacc::AcceleratorBuffer>> &, const int) const;

  mutable std::vector<ParametrizedCircuitLayer> m_layers;
  // Keeps track of the term and the index in the kernel sequence.
  mutable std::unordered_map<std::string, size_t> m_metricTermToIdx;

public:
  bool initialize(const HeterogeneousMap &parameters) override;
  const std::vector<std::string> requiredParameters() const override;
  void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
  const std::string name() const override { return "pds-vqs"; }
  const std::string description() const override { return ""; }

  DEFINE_ALGORITHM_CLONE(PDS_VQS)
};

class RiemannianMetricOptimizer : public Optimizer {
  OptResult optimize(OptFunction &) override;
  const bool isGradientBased() const override { return true; }
  const std::string name() const override { return "riemann"; }
  const std::string description() const override { return ""; }
};

//
} // namespace algorithm
} // namespace xacc

#endif