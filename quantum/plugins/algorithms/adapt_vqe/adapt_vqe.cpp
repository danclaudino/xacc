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
#include "xacc_observable.hpp"
#include "Circuit.hpp"

#include <memory>
#include <iomanip>

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

  if (parameters.keyExists<int>("maxiter")) {
    _maxIter = parameters.get<int>("maxiter");
  } 
  
  if (parameters.keyExists<double>("threshold")) {
    _threshold = parameters.get<double>("threshold");
  }

  if (parameters.keyExists<double>("print_threshold")) {
    _printThreshold = parameters.get<double>("print_threshold");
  }

  pool = parameters.getString("pool");

  // Check if Observable is Fermion or Pauli and manipulate accordingly
  if (observable->toString().find("^") != std::string::npos){

    auto jw = xacc::getService<ObservableTransform>("jw");
    if (std::dynamic_pointer_cast<FermionOperator>(observable)) {
      observable = jw->transform(observable);
    } else {
      auto fermionObservable = xacc::quantum::getObservable("fermion", observable->toString());
      observable = jw->transform(std::dynamic_pointer_cast<Observable>(fermionObservable));      
    }
  
  } else if (observable->toString().find("X") != std::string::npos
            || observable->toString().find("Y") != std::string::npos
            || observable->toString().find("Z") != std::string::npos
            && !std::dynamic_pointer_cast<PauliOperator>(observable)){

    // observable is PauliOperator, but does not cast down to it
    // Not sure about the likelyhood of this happening, but want to cover all bases
    auto pauliObservable = xacc::quantum::getObservable("pauli", observable->toString());
    observable = std::dynamic_pointer_cast<Observable>(pauliObservable);

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
  std::vector<std::shared_ptr<Observable>> pauliOps, ansatzOperators;
  auto jw = xacc::getService<ObservableTransform>("jw");
  

  // Mean-field state
  std::size_t j;
  for (int i = 0; i < nElectrons/2; i++) {
    j = (std::size_t) i;
    auto alphaXGate = ansatzRegistry->createInstruction("X", std::vector<std::size_t>{j});
    ansatzInstructions->addInstruction(alphaXGate);
    j = (std::size_t) (i + buffer->size()/2);
    auto betaXGate = ansatzRegistry->createInstruction("X", std::vector<std::size_t>{j});
    ansatzInstructions->addInstruction(betaXGate);
  }

  // Vector of non-vanishing commutators, need to compute them only once
  std::vector<std::shared_ptr<Observable>> commutators;
  for (auto op : operators){
    std::shared_ptr<Observable> pauliOp;
    if(std::dynamic_pointer_cast<PauliOperator>(op)){
      pauliOp = op;
    } else {
      pauliOp = jw->transform(op);
    }
    if (std::dynamic_pointer_cast<PauliOperator>(pauliOp)->getTerms().size() != 0){
      pauliOps.push_back(pauliOp);
      commutators.push_back(observable->commutator(pauliOp));
    }
  }

  std::vector<double> x; // these are the variational parameters
  double oldEnergy = 0.0;

  std::cout << "Operator pool: " << operatorPool->name() << "\n";
  std::cout << "Number of operators in the pool: " << operators.size() << "\n\n";

  for (int iter = 0; iter < _maxIter; iter++){

    std::cout << "Iteration: " << iter + 1 << "\n";
    std::cout << "Computing [H, A]\n" << "\n";
    std::cout << "Printing commutators with absolute value above " << _printThreshold << "\n";

    int maxCommutatorIdx = 0;
    double maxCommutator = 0.0;
    double gradientNorm = 0.0;

    // Loop over non-vanishing commutators
    for (int operatorIdx = 0; operatorIdx < commutators.size(); operatorIdx++){

      // observe the commutators with the updated circuit ansatz
      auto grad_vqe = xacc::getAlgorithm(
          "vqe", {std::make_pair("observable", commutators[operatorIdx]),
                  std::make_pair("optimizer", optimizer),
                  std::make_pair("accelerator", accelerator),
                  std::make_pair("ansatz", ansatzInstructions)});
      auto tmp_buffer = xacc::qalloc(buffer->size());
      auto commutatorValue = grad_vqe->execute(tmp_buffer, x)[0];

      if(abs(commutatorValue) > _printThreshold){
        std::cout << std::setprecision(12) << "[H," << operatorIdx << "] = " << commutatorValue << "\n";
      }

      if(abs(commutatorValue) > maxCommutator){
        maxCommutatorIdx = operatorIdx;
        maxCommutator = abs(commutatorValue);
      }

      gradientNorm += commutatorValue * commutatorValue;
    }

    gradientNorm = std::sqrt(gradientNorm);
    std::cout << "Max gradient component: [H, " << maxCommutatorIdx << "] = " << maxCommutator << " a.u.\n";
    std::cout << "Norm of gradient vector: " << gradientNorm << " a.u.\n";

    if (gradientNorm < _threshold) {
      std::cout << "\nADAPT-VQE converged in " << iter << " iterations.\n";
      std::cout << "ADAPT-VQE energy: " << oldEnergy << " a.u.\n";
      std::cout << "Optimal parameters: \n";
      for (auto param : x){
        std::cout << param << "\n";
      }
      std::cout << "\n";
      return; 

    } else if (iter < _maxIter) {

      ansatzOperators.push_back(pauliOps[maxCommutatorIdx]);

      // Instruction service for the operator to be added to the ansatz
      auto maxCommutatorGate = std::dynamic_pointer_cast<quantum::Circuit>(
          xacc::getService<Instruction>("exp_i_theta"));

      // Create instruction for new operator
      maxCommutatorGate->expand(
          {std::make_pair("pauli", pauliOps[maxCommutatorIdx]->toString()),
          std::make_pair("param_id", std::string("x") + std::to_string(iter)),
          std::make_pair("no-i", true)});

      // Add instruction for new operator to the current ansatz
      ansatzInstructions->addVariable(std::string("x") + std::to_string(iter));
      for (auto& inst : maxCommutatorGate->getInstructions()){  
        ansatzInstructions->addInstruction(inst);
      }

      // Call VQE optimization 
      auto sub_vqe = xacc::getAlgorithm(
          "vqe", {std::make_pair("observable", observable),
                  std::make_pair("optimizer", optimizer),
                  std::make_pair("accelerator", accelerator),
                  std::make_pair("ansatz", ansatzInstructions)});
      sub_vqe->execute(buffer);

      auto newEnergy = (*buffer)["opt-val"].as<double>();
      x = (*buffer)["opt-params"].as<std::vector<double>>();
      oldEnergy = newEnergy;

      std::cout << "\nEnergy at iteration " << iter + 1 << ": " << newEnergy << "\n";
      std::cout << "Parameters at iteration " << iter + 1 << ": \n";
      for (auto param : x){
        std::cout << param << "\n";
      }
      std::cout << "\n";

    } else {
      std::cout << "ADAPT-VQE did not converge in " << _maxIter << " iterations.\n";
      return;
    }
      
  }
}

} // namespace adapt_vqe
} // namespace xacc