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
#ifndef XACC_ADAPT_VQE_GRADIENT_HPP_
#define XACC_ADAPT_VQE_GRADIENT_HPP_

#include "adapt_vqe.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "AlgorithmGradientStrategy.hpp"
#include <string>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {
namespace algorithm {

class ADAPT_VQE_Gradient : public AlgorithmGradientStrategy {

protected:

  std::vector<std::shared_ptr<Observable>> ops; // operator in the current ansatz
  std::shared_ptr<Observable> H; // Hamiltonian observable
  std::vector<int> nInstructionsElement; // # of instructions for each element in gradient vector
  std::vector<double> coefficients; // coefficient that multiplies Pauli term

public:

  bool optionalParameters(const HeterogeneousMap parameters) override {

    if (!parameters.keyExists<std::shared_ptr<Observable>>("observable")){
      std::cout << "ADAPT-VQE gradient computation requires observable.\n"; 
      return false;
    }

    if (!parameters.keyExists<std::vector<std::shared_ptr<Observable>>>("operators")){
      std::cout << "ADAPT-VQE gradient computation requires excitation operators.\n"; 
      return false;
    }

    H = parameters.get<std::shared_ptr<Observable>>("observable");
    ops = parameters.get<std::vector<std::shared_ptr<Observable>>>("operators");

    return true;
    
  }

  std::vector<std::shared_ptr<CompositeInstruction>>
  getGradientExecutions(std::shared_ptr<CompositeInstruction> circuit, const std::vector<double> &x) override {

    // lambda to implement the commutator expansion up to 2nd order
    //std::complex<double> i(0, 1.0);
    auto commutator2ndOrder = [&](const std::vector<double> x, const int idx) {

      
      auto Ti = std::dynamic_pointer_cast<PauliOperator>(ops[idx]);
      for (int k = idx + 1; k < x.size(); k++){
        auto comm1stOrder = std::dynamic_pointer_cast<PauliOperator>(ops[k]->commutator(Ti));
        //std::cout << comm1stOrder->toString() <<"\n";
        comm1stOrder->operator*=(-0.5 * x[k]);
        Ti->operator+=(*comm1stOrder);
        
        auto comm2stOrder = std::dynamic_pointer_cast<PauliOperator>(ops[k]->commutator(comm1stOrder));
        comm2stOrder->operator*=(x[k] * x[k] / 8.0);
        Ti->operator+=(*comm2stOrder);

      }
      return Ti;                                  
    };

    // the derivative of latest added operator is simply [H, A]
    std::vector<std::shared_ptr<Observable>> gradientVector;
    gradientVector.push_back(H->commutator(ops.back())); 
    
    // for the remaining operators
    // [H, e^N e^(N-1)...e^(k+1) Tk e^-(k+1)...e^-(N-1) e^(-N)]
    // TODO fix this
    for (int i = 0; i < ops.size() - 1; i++){

      auto nestedCommutators = std::dynamic_pointer_cast<Observable>(commutator2ndOrder(x, i));
      auto grad = H->commutator(nestedCommutators);
      gradientVector.insert(gradientVector.begin(), std::dynamic_pointer_cast<Observable>(grad));

    }

    // Now loop over the gradient operators, observe them and store the instructions
    std::vector<std::shared_ptr<CompositeInstruction>> gradientInstructions;

    for(auto g : x){std::cout << "x " << g <<"\n";}
    std::vector<int> tmp(ops.size());
    for (int i = 0; i < gradientVector.size(); i++){

      auto kernels = gradientVector[i]->observe(circuit);

      double identityCoeff = 0.0;
      for (auto &f : kernels) {
        std::complex<double> coeff = f->getCoefficient();

        int nFunctionInstructions = 0;
        if (f->getInstruction(0)->isComposite()) {
          nFunctionInstructions =
              circuit->nInstructions() + f->nInstructions() - 1;
        } else {
          nFunctionInstructions = f->nInstructions();
        }

        if (nFunctionInstructions > circuit->nInstructions()) {
          auto evaled = f->operator()(x);
          gradientInstructions.push_back(evaled);
          coefficients.push_back(std::real(coeff));
        } else {
          identityCoeff += std::real(coeff);
        }

      }

      tmp[i] = gradientInstructions.size();

    }

    nInstructionsElement = tmp;
    return gradientInstructions;

  }

  void compute(std::vector<double> &dx, std::vector<std::shared_ptr<AcceleratorBuffer>> results) override {

    for (int gradTerm = 0; gradTerm < dx.size(); gradTerm++){ // loops over the number of entries in the gradient vector

      double gradElement = 0.0;
      for (int instElement = (gradTerm == 0) ? 0 : nInstructionsElement[gradTerm - 1];
          instElement < nInstructionsElement[gradTerm];
          instElement++) { //compute value of derivative for gradTerm entry
        auto expval = results[instElement]->getExpectationValueZ();
        gradElement += expval * coefficients[instElement];
        //if(gradTerm == 1){
        //std::cout << instElement << " exp " << expval << "coeff " << coefficients[instElement] <<"\n";
        //}
      }
      dx[gradTerm] = -gradElement / 2.0;
    }
    coefficients.clear();

    for(auto g : dx){std::cout << "grad " << g <<"\n";}

    std::reverse(dx.begin(), dx.end());
    return;
  }

  const std::string name() const override { return "adapt-vqe-gradient"; }
  const std::string description() const override { return ""; }

};

}
}

#endif