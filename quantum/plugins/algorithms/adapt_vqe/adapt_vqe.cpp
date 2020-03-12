/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
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
#include "xacc.hpp"
#include "xacc_service.hpp"

#include <memory>
#include <iomanip>
#include "<cmath>"

using namespace xacc;

namespace xacc {
namespace algorithm {
    
bool ADAPT_VQE::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    std::cout << "Obs was false\n";
    return false;
  } else if (!parameters.pointerLikeExists<CompositeInstruction>(
                 "ansatz")) {
    std::cout << "Ansatz was false\n";
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
  kernel = parameters.get<std::shared_ptr<CompositeInstruction>>("ansatz");
  nElectrons = parameters.get<int>("nElectrons");

  if (parameters.stringExists("maxiter")) {
    _maxIter = parameters.get<int>("maxiter");
  }
  //observable = parameters.getPointerLike<Observable>("observable");
  //optimizer = parameters.getPointerLike<Optimizer>("optimizer");
  //accelerator = parameters.getPointerLike<Accelerator>("accelerator");
  //kernel = parameters.getPointerLike<CompositeInstruction>("ansatz");
  pool = parameters.getString("pool");
  
  return true;
}

const std::vector<std::string> ADAPT_VQE::requiredParameters() const {
  return {"observable", "optimizer", "accelerator", "ansatz", "nElectrons", "pool"};
}

void ADAPT_VQE::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  // Here we just need to make a lambda kernel
  // to optimize that makes calls to the targeted QPU.
  
  OptFunction f(
      [&, this](const std::vector<double> &x, std::vector<double> &dx) {
        std::vector<double> coefficients;
        std::vector<std::string> kernelNames;
        std::vector<std::shared_ptr<CompositeInstruction>> fsToExec;

        double identityCoeff = 0.0;
        for (auto &f : kernels) {
          kernelNames.push_back(f->name());
          std::complex<double> coeff = f->getCoefficient();

          int nFunctionInstructions = 0;
          if (f->getInstruction(0)->isComposite()) {
            nFunctionInstructions =
                kernel->nInstructions() + f->nInstructions() - 1;
          } else {
            nFunctionInstructions = f->nInstructions();
          }

          if (nFunctionInstructions > kernel->nInstructions()) {
            auto evaled = f->operator()(x);
            fsToExec.push_back(evaled);
            coefficients.push_back(std::real(coeff));
          } else {
            identityCoeff += std::real(coeff);
          }
        }

        auto tmpBuffer = xacc::qalloc(buffer->size());
        accelerator->execute(tmpBuffer, fsToExec);
        auto buffers = tmpBuffer->getChildren();

        double energy = identityCoeff;
        auto idBuffer = xacc::qalloc(buffer->size());
        idBuffer->addExtraInfo("coefficient", identityCoeff);
        idBuffer->setName("I");
        idBuffer->addExtraInfo("kernel", "I");
        idBuffer->addExtraInfo("parameters", x);
        idBuffer->addExtraInfo("exp-val-z", 1.0);
        if (accelerator->name() == "ro-error")
            idBuffer->addExtraInfo("ro-fixed-exp-val-z", 1.0);
        buffer->appendChild("I", idBuffer);

        if (buffers[0]->hasExtraInfoKey(
                "purified-energy")) { // FIXME Hack for now...
          energy = buffers[0]->getInformation("purified-energy").as<double>();
          for (auto &b : buffers) {
            // b->addExtraInfo("parameters", initial_params);
            buffer->appendChild(b->name(), b);
          }
        } else {
          for (int i = 0; i < buffers.size(); i++) {
            auto expval = buffers[i]->getExpectationValueZ();
            energy += expval * coefficients[i];
            buffers[i]->addExtraInfo("coefficient", coefficients[i]);
            buffers[i]->addExtraInfo("kernel", fsToExec[i]->name());
            buffers[i]->addExtraInfo("exp-val-z", expval);
            buffers[i]->addExtraInfo("parameters", x);
            buffer->appendChild(fsToExec[i]->name(), buffers[i]);
          }
        }

        std::stringstream ss;
        ss << "E(" << ( !x.empty() ? std::to_string(x[0]) : "");
        for (int i = 1; i < x.size(); i++)
          ss << "," << x[i];
        ss << ") = " << std::setprecision(12) << energy;
        xacc::info(ss.str());
        return result;
      },
      kernel->nVariables());

  std::shared_ptr<CompositeInstruction>
  circuit(const std::unordered_map<std::string, Term>){

    auto gateRegistry = xacc::getIRProvider("quantum");
    for (auto inst : terms) {

      // for (auto inst : s) {
      Term spinInst = inst.second;

      // Get the individual pauli terms
      auto termsMap = std::get<2>(spinInst);

      std::vector<std::pair<int, std::string>> terms;
      for (auto &kv : termsMap) {
        if (kv.second != "I" && !kv.second.empty()) {
          terms.push_back({kv.first, kv.second});
        }
      }
      // The largest qubit index is on the last term
      int largestQbitIdx = terms[terms.size() - 1].first;
      auto tempFunction = gateRegistry->createComposite("temp", {spinInst.var()});

      for (int i = 0; i < terms.size(); i++) {

        std::size_t qbitIdx = terms[i].first;
        auto gateName = terms[i].second;

        if (i < terms.size() - 1) {
          std::size_t tmp = terms[i + 1].first;
          auto cnot = gateRegistry->createInstruction(
              "CNOT", std::vector<std::size_t>{qbitIdx, tmp});
          tempFunction->addInstruction(cnot);
        }

        if (gateName == "X") {
          auto hadamard = gateRegistry->createInstruction(
              "H", std::vector<std::size_t>{qbitIdx});
          tempFunction->insertInstruction(0, hadamard);
        } else if (gateName == "Y") {
          auto rx = gateRegistry->createInstruction(
              "Rx", std::vector<std::size_t>{qbitIdx});
          InstructionParameter p(pi / 2.0);
          rx->setParameter(0, p);
          tempFunction->insertInstruction(0, rx);
        }

        // Add the Rotation for the last term
        if (i == terms.size() - 1) {
          // FIXME DONT FORGET DIVIDE BY 2
          std::stringstream ss;
          ss << 2 * std::imag(std::get<0>(spinInst)) << " * "
            << std::get<1>(spinInst);
          auto rz = gateRegistry->createInstruction(
              "Rz", std::vector<std::size_t>{qbitIdx});

          InstructionParameter p(ss.str());
          rz->setParameter(0, p);
          tempFunction->addInstruction(rz);
        }
      }

      int counter = tempFunction->nInstructions();
      // Add the instruction on the backend of the circuit
      for (int i = terms.size() - 1; i >= 0; i--) {

        std::size_t qbitIdx = terms[i].first;
        auto gateName = terms[i].second;

        if (i < terms.size() - 1) {
          std::size_t tmp = terms[i + 1].first;
          auto cnot = gateRegistry->createInstruction(
              "CNOT", std::vector<std::size_t>{qbitIdx, tmp});
          tempFunction->insertInstruction(counter, cnot);
          counter++;
        }

        if (gateName == "X") {
          auto hadamard = gateRegistry->createInstruction(
              "H", std::vector<std::size_t>{qbitIdx});
          tempFunction->addInstruction(hadamard);
        } else if (gateName == "Y") {
          auto rx = gateRegistry->createInstruction(
              "Rx", std::vector<std::size_t>{qbitIdx});
          InstructionParameter p(-1.0 * (pi / 2.0));
          rx->setParameter(0, p);
          tempFunction->addInstruction(rx);
        }
      }
      // Add to the total UCCSD State Prep function
      for (auto inst : tempFunction->getInstructions()) {
        addInstruction(inst);
      }
    }
/*
    for (int i = nElectrons - 1; i >= 0; i--) {
      std::size_t j = (std::size_t)i;
      auto xGate =
          gateRegistry->createInstruction("X", std::vector<std::size_t>{j});
      insertInstruction(0, xGate);
    }
*/
    return gateRegistry;
  }

  // auto kernels = observable->observe(xacc::as_shared_ptr(kernel));
  // create circuit for HF/MF state assuming all ups then all downs
  //auto provider = xacc::getIRProvider("quantum");

  auto kernel = provider->createComposite("circuit");

  for (size_t i; i < nElectrons/2; i++){

    auto xUp = provider->createInstruction("X", {i});
    auto xDown = provider->createInstruction("X", {i + nElectrons});
    kernel->addInstructions({xUp, xDown});

  }

  auto operatorPool = xacc::getService<OperatorPool>(pool);
  auto operators = operatorPool->generate(buffer->size(), nElectrons);
  std::unordered_map<std::string, Term> pauliOperators =
    std::dynamic_pointer_cast<PauliOperator>(operators)->getTerms();
  std::unordered_map<std::string, Term> ansatzTerms;

  std::vector<double> x; // these are the thetas
  // auto thetaIterator = theta.insert(theta.begin(), 0.0);
  for (int iter; iter < _maxIter; iter++){

    // compute commutators
    int largestCommutatorIdx = -1;
    double largestCommutator = 0.0;
    double oldEnergy = 0.0;
    for (int operatorIdx; operatorIdx < terms.size(); operatorIdx++){

      //compute commutator for operatorIdx
      auto operatorCommutator = commutator();
      if(operatorCommutator < largestCommutator){
        largestCommutatorIdx = operatorIdx;
        largestCommutator = operatorCommutator;
      }
    }

    // add operator with largest commutator to the ansatz
    ansatzTerms.push_back(pauliOperators[largestCommutatorIdx]);
    auto ansatzCircuit = circuit(ansatzTerms);
    kernel->addInstruction(ansatzCircuit)
    // call optimize
    auto newResult = optimizer->optimize(f);
    auto newEnergy = newResult.first;
    if (abs(newEnergy - oldEnergy) <= threshold){
      std::cout << "ADAPT-VQE converged in" << _maxIter << "iterations.\n";
      return true;
    } else if (abs(newEnergy - oldEnergy) > threshold && iter < _maxIter){
      oldEnergy = newEnergy;
      x.push_back(0.0);
    } else {
      std::cout << "ADAPT-VQE did not converge in" << _maxIter << " iterations.\n";
      return false;
    }
      
  }

  buffer->addExtraInfo("opt-val", ExtraInfo(result.first));
  buffer->addExtraInfo("opt-params", ExtraInfo(result.second));
  return;
}
} // namespace adapt_vqe
} // namespace xacc
