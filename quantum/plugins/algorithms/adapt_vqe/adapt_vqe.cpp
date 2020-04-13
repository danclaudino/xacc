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

  // Check if Observable is Fermion or Pauli and manipulate the pointer accordingly
if (!std::dynamic_pointer_cast<FermionOperator>(observable) && observable->toString().find("^") != std::string::npos) {
    // Case 1: observable is FermionOperator, but does not cast down to it
    std::cout << observable->toString();
    auto fermionObservable = xacc::quantum::getObservable("fermion", std::string("(-0.165494,-0)  2^ 1^ 2 1 "));
    //observable = jw->transform(fermionObservable);
    observable = std::dynamic_pointer_cast<Observable>(fermionObservable);
  
  } else if (observable->toString().find("X") != std::string::npos
            || observable->toString().find("Y") != std::string::npos
            || observable->toString().find("Z") != std::string::npos
            && !std::dynamic_pointer_cast<PauliOperator>(observable)){
    // Case 2: observable is PauliOperator, but does not cast down to it
    // Not sure about the likelyhood of this happening, but want to cover all bases
    auto pauliObservable = xacc::quantum::getObservable("pauli", observable->toString());
    observable = std::dynamic_pointer_cast<Observable>(pauliObservable);

  } // if Obser casts down, nothing is neededg


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
  
  // instructions for mean-field state
  for (int i = nElectrons - 1; i >= 0; i--) {
    std::size_t j = (std::size_t)i;
    auto xGate =
      ansatzRegistry->createInstruction("X", std::vector<std::size_t>{j});
    ansatzInstructions->addInstruction(xGate);
  }

//std::cout << observable.asPauli->toString();
  std::vector<double> x; // these are the thetas
  double oldEnergy = 0.0;
  std::cout << "Observable\n\n" << observable->toString();
        if(std::dynamic_pointer_cast<FermionOperator>(observable)){
        std::cout << "Observable is fermion\n\n";
      } else if (std::dynamic_pointer_cast<PauliOperator>(observable)) {
        std::cout << "Observable is Pauli\n\n";
      } else {
std::cout << "Observable is Observable \n\n";
      }
  for (int iter = 0; iter < 2; iter++){


    std::cout << "Iteration: " << iter + 1 << std::endl;
    std::cout << "Computing [H, A]\n" << std::endl;

    // compute commutators
    int largestCommutatorIdx = 0;
    double largestCommutator = 0.0;
    double gradientNorm = 0.0;;
    
    for (int operatorIdx = 0; operatorIdx < operators.size(); operatorIdx++){

      //compute commutator for operatorIdx
      //PauliOperator& H = *std::dynamic_pointer_cast<PauliOperator>(obs).get();
      //PauliOperator& A = *std::dynamic_pointer_cast<PauliOperator>(operators[operatorIdx]).get();
      //PauliOperator commutator = H * A - A * H;
      //auto operatorCommutatorPtr = std::shared_ptr<Observable>(&commutator, [](Observable *) {});


      
      auto operatorCommutatorPtr = observable->commutator(operators[operatorIdx]);
      
      auto grad_vqe_energy = xacc::getAlgorithm(
          "vqe", {std::make_pair("observable", operatorCommutatorPtr),
                  std::make_pair("optimizer", optimizer),
                  std::make_pair("accelerator", accelerator),
                  std::make_pair("ansatz", ansatzInstructions)});

      auto tmp_buffer = xacc::qalloc(buffer->size());
      auto commutatorValue = grad_vqe_energy->execute(tmp_buffer, x)[0];
      std::cout << "[H," << operatorIdx << "] = " << commutatorValue << std::endl;
    
      if(abs(commutatorValue) > largestCommutator){
        largestCommutatorIdx = operatorIdx;
        largestCommutator = commutatorValue;
      }

      gradientNorm += commutatorValue * commutatorValue;
    }

    gradientNorm = std::sqrt(gradientNorm);
    std::cout << "Norm of gradient vector: " << gradientNorm << "a.u.\n";

    if (gradientNorm < _threshold) {
      std::cout << "ADAPT-VQE converged in " << iter << " iterations.\n";
      std::cout << "ADAPT-VQE energy:" << oldEnergy << " a.u.\n";
      return; 

    } else if (iter < _maxIter) {
      auto largestCommutatorGate = std::dynamic_pointer_cast<quantum::Circuit>(
          xacc::getService<Instruction>("exp_i_theta"));

      largestCommutatorGate->expand(
          {std::make_pair("pauli", operators[largestCommutatorIdx]->toString()),
          std::make_pair("param_id", std::string("x") + std::to_string(iter)),
          std::make_pair("no-i", true)});

      //auto adaptAnsatz = ansatzInstructions->clone();
      ansatzInstructions->addVariable(std::string("x") + std::to_string(iter));
      for (auto& inst : largestCommutatorGate->getInstructions()){  
        ansatzInstructions->addInstruction(inst);
      }

      auto sub_vqe = xacc::getAlgorithm(
          "vqe", {std::make_pair("observable", observable),
                  std::make_pair("optimizer", optimizer),
                  std::make_pair("accelerator", accelerator),
                  std::make_pair("ansatz", ansatzInstructions)});
      sub_vqe->execute(buffer);

      auto newEnergy = (*buffer)["opt-val"].as<double>();
      auto newParams = (*buffer)["opt-params"].as<std::vector<double>>();

      for (int i = 0; i < newParams.size() - 1; i++){
        x.push_back(newParams[i]);
      }
      x.push_back(newParams.back());

      std::cout << "Energy at iteration " << iter + 1 << " : " << newEnergy << "\n\n" << std::endl;
      oldEnergy = newEnergy;

    } else {
      std::cout << "ADAPT-VQE did not converge in" << _maxIter << " iterations.\n";
      return;
    }
      
  }
}

} // namespace adapt_vqe
} // namespace xacc