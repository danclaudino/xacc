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

#include "adapt_vqe.hpp"

#include "Observable.hpp"
#include "PauliOperator.hpp"
#include "FermionOperator.hpp"
#include "ObservableTransform.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Circuit.hpp"

#include <memory>
#include <iomanip>
#include <cmath>
#include <typeinfo>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {
namespace algorithm {
    
bool ADAPT_VQE::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    std::cout << "Obs was false\n";
    return false;
  } else if (!parameters.pointerLikeExists<Accelerator>(
                 "accelerator")) {
    std::cout << "Acc was false\n";
    return false;
  } else if(!parameters.stringExists("pool")){
    return false;
  }

  optimizer = parameters.get<std::shared_ptr<Optimizer>>("optimizer");
  observable = parameters.get<std::shared_ptr<Observable>>("observable");
  accelerator = parameters.get<std::shared_ptr<Accelerator>>("accelerator");
  nElectrons = parameters.get<int>("nElectrons");

  if (parameters.stringExists("maxiter")) {
    _maxIter = parameters.get<int>("maxiter");
  } 
  
  if (parameters.stringExists("threshold")) {
    _threshold = parameters.get<double>("threshold");
  } 
  pool = parameters.getString("pool");
  
  if (std::dynamic_pointer_cast<FermionOperator>(observable)) {
    auto jw = xacc::getService<ObservableTransform>("jw");
    observable = jw->transform(observable);
  } 

  return true;
}

const std::vector<std::string> ADAPT_VQE::requiredParameters() const {
  return {"observable", "optimizer", "accelerator", "nElectrons", "pool"};
}

void ADAPT_VQE::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  auto ansatzRegistry = xacc::getIRProvider("quantum");
  auto ansatzInstructions = ansatzRegistry->createComposite("ansatzCircuit");
  auto operatorPool = xacc::getService<OperatorPool>(pool);
  auto operators = operatorPool->generate(buffer->size(), nElectrons);
  auto vqe = xacc::getAlgorithm("vqe");
  auto largestCommutatorGate = std::dynamic_pointer_cast<quantum::Circuit>(
    xacc::getService<Instruction>("exp_i_theta"));

  HeterogeneousMap vqeParameters;
  vqeParameters.insert("optimizer", optimizer);
  vqeParameters.insert("accelerator", accelerator);
  
  // instructions for mean-field state
  for (int i = nElectrons - 1; i >= 0; i--) {
    std::size_t j = (std::size_t)i;
    auto xGate =
      ansatzRegistry->createInstruction("X", std::vector<std::size_t>{j});
    ansatzInstructions->addInstruction(xGate);
  }

  std::vector<double> x; // these are the thetas
  //x.push_back(0.0);
  double oldEnergy = 0.0;
  for (int iter = 0; iter < 1; iter++){

    std::cout << "Iteration: " << iter << std::endl;
    std::cout << "Computing [H, A]" << std::endl;

    // compute commutators
    int largestCommutatorIdx = 0;
    double largestCommutator = 0.0;
    for (int operatorIdx = 0; operatorIdx < operators.size(); operatorIdx++){

      //compute commutator for operatorIdx
      PauliOperator& H = *std::dynamic_pointer_cast<PauliOperator>(observable).get();
      PauliOperator& A = *std::dynamic_pointer_cast<PauliOperator>(operators[operatorIdx]).get();
      PauliOperator commutator = H * A - A * H;
      auto operatorCommutatorPtr = std::shared_ptr<Observable>(&commutator, [](Observable *) {});

      vqeParameters.insert("observable", operatorCommutatorPtr);
      vqeParameters.insert("ansatz", ansatzInstructions);
      vqe->initialize(vqeParameters);

      auto commutatorValue = vqe->execute(buffer, x)[0];
      std::cout << "[H," << operatorIdx << "] = " << commutatorValue << std::endl;
    
      if(commutatorValue > largestCommutator){
        largestCommutatorIdx = operatorIdx;
        largestCommutator = commutatorValue;
      } 
    }

    largestCommutatorGate->expand({std::make_pair("pauli", operators[largestCommutatorIdx]->toString()), std::make_pair("param_id", std::string("x")+std::to_string(iter))});
    ansatzInstructions->addVariable(std::string("x") + std::to_string(iter));
    for (auto& inst : largestCommutatorGate->getInstructions()){  
      ansatzInstructions->addInstruction(inst);
    }
    vqeParameters.insert("observable", observable);
    vqeParameters.insert("ansatz", ansatzInstructions);
    vqe->initialize(vqeParameters);
    vqe->execute(buffer);

    auto newEnergy = (*buffer)["opt-val"].as<double>();
    std::cout << "Energy at iteration " << iter << " : " << newEnergy << std::endl;
    if (abs(newEnergy - oldEnergy) <= _threshold){
      std::cout << "ADAPT-VQE converged in " << iter << " iterations.\n";
      std::cout << "ADAPT-VQE energy:" << newEnergy << "a.u.\n";
      return;
    } else if (abs(newEnergy - oldEnergy) > _threshold && iter < _maxIter){
      oldEnergy = newEnergy;
      x.push_back(0.0);
    } else {
      std::cout << "ADAPT-VQE did not converge in" << _maxIter << " iterations.\n";
      return;
    }
      
  }
}

} // namespace adapt_vqe
} // namespace xacc