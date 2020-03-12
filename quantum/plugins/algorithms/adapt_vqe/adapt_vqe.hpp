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
// include "Observable.hpp"
#include <vector>

namespace xacc{
namespace algorithm{

class OperatorPool : public Identifiable {
public:
  virtual bool isValidOperatorPool(const std::string &operatorPool) = 0;

  virtual std::string 
  constructOperatorString(const int p, const int q, const std::string norm = "1.0");
  virtual std::string
  constructOperatorString(const int p, const int q, const int r, 
    const int s, const std::string norm = "1.0");

  virtual std::vector<std::shared_ptr<Observable>>
  generate(const int &nQubits, const int &nElectrons);

};

class ADAPT_VQE : public Algorithm {
protected:
  std::shared_ptr<Observable> observable;
  std::shared_ptr<Optimizer> optimizer;
  std::shared_ptr<CompositeInstruction> kernel; // keep this for now
  std::shared_ptr<Accelerator> accelerator;
  int nElectrons;
  std::string pool;
  int _maxIter;
  // std::vector<double> initial_params;
  // std::string gradient; will eventually come back here

  HeterogeneousMap _parameters;

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