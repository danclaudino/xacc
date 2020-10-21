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
#ifndef XACC_ALGORITHM_RBM_CHEMISTRY_HPP_
#define XACC_ALGORITHM_RBM_CHEMISTRY_HPP_

#include "Algorithm.hpp"
#include <Eigen/Dense>
#include <memory>
#include <vector>
#include "AnnealingProgram.hpp"
#include "EmbeddingAlgorithm.hpp"
#include "CompositeInstruction.hpp"

using bitStringPair = std::pair<std::string, std::string>;
using Counts = std::map<std::string, int>;

namespace xacc {
namespace algorithm {
class RBMChemistry : public Algorithm {

protected:
  // std::shared_ptr<xacc::quantum::AnnealingProgram> rbm;
  std::shared_ptr<CompositeInstruction> rbm;
  Observable *observable;
  Accelerator *accelerator;
  HeterogeneousMap _parameters;
  int nHidden, nVisible;
  std::string default_emb_algo = "cmr";

  std::vector<std::pair<bitStringPair, double>> hamiltonian;

  double wInitMax = 0.1, eps = 1.0e-4, minCoefficient = 1.0e-4;

  std::vector<double> random_vector(const double l_range, const double r_range,
                                    const std::size_t size) const;

  std::string generateBitString(unsigned decimal) const;

  Eigen::VectorXd computeUpdate(std::vector<double> wb, Counts counts) const;

  double computeScale(const std::vector<double> a, const std::vector<double> b,
                      const std::vector<double> w,
                      const std::map<int, std::vector<int>> embedding) const;

  Counts sample(std::shared_ptr<AcceleratorBuffer> &buffer, const std::shared_ptr<CompositeInstruction> problem) const;

public:
  bool initialize(const HeterogeneousMap &parameters) override;
  const std::vector<std::string> requiredParameters() const override;

  void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
  const std::string name() const override { return "rbm-chemistry"; }
  const std::string description() const override { return ""; }
  DEFINE_ALGORITHM_CLONE(RBMChemistry)
};
} // namespace algorithm
} // namespace xacc
#endif
