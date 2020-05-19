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
#ifndef XACC_ALGORITHM_ADAPT_VQE_HPP_
#define XACC_ALGORITHM_ADAPT_VQE_HPP_

#include "Algorithm.hpp"
#include "Observable.hpp"
#include "PauliOperator.hpp"

using namespace xacc::quantum;

namespace xacc{
namespace algorithm{

class OperatorPool : public Identifiable {

public:
  virtual bool isValidOperatorPool(const std::string &operatorPool) = 0;

  virtual std::vector<std::shared_ptr<Observable>>
  generate(const int &nQubits, const int &nElectrons) = 0;

};


class GradientStrategy : public Identifiable {

  protected:
    //virtual std::vector<double> gradientElements(const std::vector<double> &x) = 0;
  virtual std::vector<std::shared_ptr<CompositeInstruction>>
  getCircuitExecutions(std::vector<std::shared_ptr<CompositeInstruction>> ansatzInstructions,
                      const std::vector<std::shared_ptr<Observable>> ansatzOperators,
                      const std::vector<double> &x) = 0;

  virtual void compute(std::vector<double> &grad, std::vector<std::shared_ptr<AcceleratorBuffer>> results) = 0;

  const std::string name() const override { return "adapt-vqe-gradient"; }
  const std::string description() const override { return ""; }

};

class ADAPT_VQE : public Algorithm {
protected:
  std::shared_ptr<Observable> observable;
  std::shared_ptr<Optimizer> optimizer;
  std::shared_ptr<Accelerator> accelerator;
  int nElectrons;
  std::string pool;
  int _maxIter = 50;
  double _threshold = 1.0e-2;
  double _printThreshold = 1.0e-10;
  std::vector<std::shared_ptr<Observable>> ansatzOperators;

  HeterogeneousMap _parameters;

  //std::vector<std::shared_ptr<Observable>> 
  //gradientVector(const std::vector<double> &x, 
              //  const std::vector<std::shared_ptr<Observable>> ops);

public:

  bool initialize(const HeterogeneousMap &parameters) override;
  const std::vector<std::string> requiredParameters() const override;
  void execute(const std::shared_ptr<AcceleratorBuffer> buffer) const override;
  const std::string name() const override { return "adapt-vqe"; }
  const std::string description() const override { return ""; }

  DEFINE_ALGORITHM_CLONE(ADAPT_VQE)

};


} // namespace algorithm
} // namespace xacc

#endif