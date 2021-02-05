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
#ifndef XACC_ALGORITHM_ADAPT_HPP_
#define XACC_ALGORITHM_ADAPT_HPP_

#include "Algorithm.hpp"
#include "PauliOperator.hpp"
#include "OperatorPool.hpp"
#include <memory>
#include <vector>

using namespace xacc::quantum;

namespace xacc {
namespace algorithm {

class ADAPT : public Algorithm {

protected:
  std::shared_ptr<Observable> observable;
  std::shared_ptr<Optimizer> optimizer;
  std::shared_ptr<Accelerator> accelerator;
  std::shared_ptr<OperatorPool> pool;
  std::shared_ptr<CompositeInstruction> initialState;
  std::string subAlgo; // sub-algorithm, either VQE or QAOA
  mutable HeterogeneousMap _parameters;

  // ADAPT parameters
  int _maxIter = 50;                // max # of ADAPT cycles // # of QAOA layers
  double _adaptThreshold = 1.0e-2;  // gradient norm threshold
  double _printThreshold = 1.0e-10; // threshold to print commutator
  double _measurementThreshold = 0.0; // threshold to ignore measurement
  bool _printOps = false; // set to true to print operators at every iteration
  int _nElectrons;        // # of electrons, used for VQE

  // indices of operators to construct initial ansatz
  std::vector<int> initialOps;
  // initial parameters for initial ansatz
  std::vector<double> initialParams;
  // name of class to compute gradient for optimization
  // defaults to parameter shift
  std::string gradStrategyName = "parameter-shift";

  // initializes new parameter from partial tomography
  // where all previous parameters are frozen
  double newParameter(const std::shared_ptr<Observable> obs,
                      const std::shared_ptr<CompositeInstruction> kernel,
                      const std::vector<double> &x, double zeroExpValue) const;

  // cache commutator measurements for a given ADAPT iteration
  mutable std::unordered_map<std::string, double> cachedMeasurements;

  // replaces vqe(buffer, {}) to enable measurement caching
  double measureOperator(const std::shared_ptr<Observable> obs,
                         const std::shared_ptr<CompositeInstruction> kernel,
                         const std::vector<double> &x,
                         const int bufferSize) const;

public:
  bool initialize(const HeterogeneousMap &parameters) override;
  const std::vector<std::string> requiredParameters() const override;
  void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
  const std::string name() const override { return "adapt"; }
  const std::string description() const override { return ""; }

  DEFINE_ALGORITHM_CLONE(ADAPT)
};

} // namespace algorithm
} // namespace xacc

#endif