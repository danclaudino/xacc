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

#include "adapt.hpp"

#include "FermionOperator.hpp"
#include "ObservableTransform.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include "Circuit.hpp"
#include "AlgorithmGradientStrategy.hpp"

#include <iomanip>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {
namespace algorithm {

bool ADAPT::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    xacc::info("Obs was false\n");
    return false;
  }

  if (!parameters.pointerLikeExists<Accelerator>("accelerator")) {
    xacc::info("Acc was false\n");
    return false;
  }

  if (!parameters.stringExists("pool")) {
    xacc::info("Pool was false\n");
    return false;
  }

  if (!parameters.stringExists("sub-algorithm")) {
    xacc::info("Sub algorithm was false\n");
    return false;
  }

  optimizer = parameters.get<std::shared_ptr<Optimizer>>("optimizer");

  observable = parameters.get<std::shared_ptr<Observable>>("observable");
  if (std::dynamic_pointer_cast<FermionOperator>(observable)) {
    observable = getService<ObservableTransform>("jw")->transform(observable);
  }

  accelerator = parameters.get<std::shared_ptr<Accelerator>>("accelerator");
  pool = xacc::getService<OperatorPool>(parameters.getString("pool"));
  subAlgo = parameters.getString("sub-algorithm");

  if (parameters.keyExists<int>("maxiter")) {
    _maxIter = parameters.get<int>("maxiter");
  }

  if (parameters.keyExists<double>("adapt-threshold")) {
    _adaptThreshold = parameters.get<double>("adapt-threshold");
  }

  if (parameters.keyExists<double>("print-threshold")) {
    _printThreshold = parameters.get<double>("print-threshold");
  }

  if (parameters.stringExists("gradient_strategy")) {
    gradStrategyName = parameters.getString("gradient_strategy");
  }

  if (parameters.pointerLikeExists<CompositeInstruction>("initial-state")) {
    initialState =
        parameters.get<std::shared_ptr<CompositeInstruction>>("initial-state");
  }

  if (parameters.keyExists<int>("n-electrons")) {
    _nElectrons = parameters.get<int>("n-electrons");
  }

  if ((subAlgo == "vqe" && !initialState) &&
      (subAlgo == "vqe" && !parameters.keyExists<int>("n-electrons"))) {

    xacc::info("VQE requires number of electrons or initial state.");
  }

  if (parameters.getString("pool") == "singlet-adapted-uccsd" &&
      parameters.keyExists<int>("n-electrons")) {

    pool->optionalParameters({{"n-electrons", _nElectrons}});
  }

  return true;
}

const std::vector<std::string> ADAPT::requiredParameters() const {
  return {"observable", "optimizer", "accelerator", "pool", "sub-algorithm"};
}

void ADAPT::execute(const std::shared_ptr<AcceleratorBuffer> buffer) const {

  std::stringstream ss;

  // Instantiate adaptive circuit and add instructions from initial state
  // initial state for chemistry VQE is usually HF
  // initial state for QAOA is all |+>
  auto provider = xacc::getIRProvider("quantum");
  auto ansatz = provider->createComposite("ansatz");

  if (initialState) {

    for (auto &inst : initialState->getInstructions()) {
      ansatz->addInstruction(inst);
    }

  } else {

    if (subAlgo == "vqe") {
      std::size_t j;
      for (int i = 0; i < _nElectrons / 2; i++) {
        j = (std::size_t)i;
        ansatz->addInstruction(provider->createInstruction("X", {j}));
        j = (std::size_t)(i + buffer->size() / 2);
        ansatz->addInstruction(provider->createInstruction("X", {j}));
      }
    }

    if (subAlgo == "QAOA") {
      std::size_t j;
      for (int i = 0; i < buffer->size(); i++) {
        j = (std::size_t)i;
        ansatz->addInstruction(provider->createInstruction("H", {j}));
      }
    }
  }

  // Generate operators in the pool
  auto operators = pool->generate(buffer->size());
  std::vector<int> ansatzOps;

  // Vector of commutators, need to compute them only once
  std::vector<std::shared_ptr<Observable>> commutators;
  if (subAlgo == "QAOA") {

    std::complex<double> i(0, 1);
    for (auto op : operators) {
      auto comm = observable->commutator(op);
      auto &tmp = *std::dynamic_pointer_cast<PauliOperator>(comm);
      tmp = tmp * i;
      auto ptr = std::dynamic_pointer_cast<Observable>(std::make_shared<PauliOperator>(tmp));
      commutators.push_back(ptr);
    }

  } else {

    for (auto op : operators) {
      commutators.push_back(observable->commutator(op));
    }

  }

  xacc::info("Operator pool: " + pool->name());
  xacc::info("Number of operators in the pool: " +
             std::to_string(operators.size()));

  double oldEnergy = 0.0;
  std::vector<double> x; // these are the variational parameters

  // start ADAPT loop
  for (int iter = 0; iter < _maxIter; iter++) {

    xacc::info("Iteration: " + std::to_string(iter + 1));
    xacc::info("Computing [H, A]");
    xacc::info("Printing commutators with absolute value above " +
               std::to_string(_printThreshold));

    if (subAlgo == "QAOA") {

      auto costHamiltonianGates = std::dynamic_pointer_cast<quantum::Circuit>(
          xacc::getService<Instruction>("exp_i_theta"));

      // Create instruction for new operator
      costHamiltonianGates->expand({{"pauli", observable->toString()},
           {"param_id", "x" + std::to_string(ansatz->nVariables())}});

      ansatz->addVariable("x" + std::to_string(ansatz->nVariables()));
      for (auto &inst : costHamiltonianGates->getInstructions()) {
        ansatz->addInstruction(inst);
      }
      x.insert(x.begin(), 0.01);
    }

    int maxCommutatorIdx = 0;
    double maxCommutator = 0.0;
    double gradientNorm = 0.0;

    // Loop over non-vanishing commutators and select the one with largest
    // magnitude
    for (int opIdx = 0; opIdx < commutators.size(); opIdx++) {

      // only compute commutators if they aren't zero
      int nTermsCommutator =
          std::dynamic_pointer_cast<PauliOperator>(commutators[opIdx])
              ->getTerms()
              .size();
      if (nTermsCommutator != 0) {

        // Print number of instructions for computing <observable>
        xacc::info("Number of instructions for commutator calculation: " +
                   std::to_string(nTermsCommutator));

        // observe the commutators with the updated circuit ansatz
        auto grad_algo = xacc::getAlgorithm(
            subAlgo, {{"observable", commutators[opIdx]},
                      {"optimizer", optimizer},
                      {"accelerator", accelerator},
                      {"ansatz", ansatz}});
        auto tmp_buffer = xacc::qalloc(buffer->size());
        auto commutatorValue = grad_algo->execute(tmp_buffer, x)[0];

        if (abs(commutatorValue) > _printThreshold) {
          ss << std::setprecision(12) << "[H," << opIdx
             << "] = " << commutatorValue;
          xacc::info(ss.str());
          ss.str(std::string());
        }

        // update maxCommutator
        if (abs(commutatorValue) > abs(maxCommutator)) {
          maxCommutatorIdx = opIdx;
          maxCommutator = commutatorValue;
        }

        gradientNorm += commutatorValue * commutatorValue;
      }
    }

    ss << std::setprecision(12) << "Max gradient component: [H, "
       << maxCommutatorIdx << "] = " << maxCommutator << " a.u.";
    xacc::info(ss.str());
    ss.str(std::string());

    gradientNorm = std::sqrt(gradientNorm);
    ss << std::setprecision(12) << "Norm of gradient vector: " << gradientNorm
       << " a.u.";
    xacc::info(ss.str());
    ss.str(std::string());

    if (gradientNorm < _adaptThreshold) { // ADAPT converged

      xacc::info("ADAPT-" + subAlgo + " converged in " + std::to_string(iter) +
                 " iterations.");
      return;

    } else if (iter < _maxIter) { // Add operator and reoptimize

      xacc::info(subAlgo + " optimization of current ansatz.");

      // keep track of growing ansatz
      ansatzOps.push_back(maxCommutatorIdx);

      // Instruction service for the operator to be added to the ansatz
      auto maxCommutatorGate = pool->getOperatorInstructions(
          maxCommutatorIdx, ansatz->nVariables());

      // Label for new variable and add it to the circuit
      ansatz->addVariable("x" + std::to_string(ansatz->nVariables()));

      // Append new instructions to current circuit
      for (auto &inst : maxCommutatorGate->getInstructions()) {
        ansatz->addInstruction(inst);
      }

      // Convergence is improved if passing initial parameters to optimizer
      // so we create a new instance of Optimizer with them
      // insert 0.0 for VQE and pi/2 for QAOA
      if (subAlgo == "vqe") {
        x.insert(x.begin(), 0.0);
      }
      if (subAlgo == "QAOA") {
        x.insert(x.begin(), xacc::constants::pi / 2.0);
      }
      auto newOptimizer = xacc::getOptimizer(optimizer->name(), 
                          {{"nlopt-optimizer", optimizer->get_algorithm()},
                           {"initial-parameters", x}});

      // Instantiate gradient class
      std::shared_ptr<AlgorithmGradientStrategy> gradientStrategy;
      if (!gradStrategyName.empty() && optimizer->isGradientBased()) {
        gradientStrategy =
            xacc::getService<AlgorithmGradientStrategy>(gradStrategyName);
        gradientStrategy->initialize({{"observable", observable}});
      } else {
        gradientStrategy = nullptr;
      }

      // Start subAlgo optimization
      auto sub_opt = xacc::getAlgorithm(
          subAlgo, {{"observable", observable},
                    {"optimizer", newOptimizer},
                    {"accelerator", accelerator},
                    {"gradient_strategy", gradientStrategy},
                    {"ansatz", ansatz}});
      sub_opt->execute(buffer);

      auto newEnergy = (*buffer)["opt-val"].as<double>();
      x = (*buffer)["opt-params"].as<std::vector<double>>();
      oldEnergy = newEnergy;

      ss << std::setprecision(12) << "Energy at ADAPT iteration " << iter + 1
         << ": " << newEnergy;
      xacc::info(ss.str());
      ss.str(std::string());

      ss << std::setprecision(12) << "Parameters at ADAPT iteration "
         << iter + 1 << ": ";
      for (auto param : x) {
        ss << param << " ";
      }
      xacc::info(ss.str());
      ss.str(std::string());

      ss << "Ansatz at ADAPT iteration " << std::to_string(iter + 1) << ": ";
      for (auto op : ansatzOps) {
        ss << op << " ";
      }
      xacc::info(ss.str());
      ss.str(std::string());

      buffer->addExtraInfo("opt-val", ExtraInfo(oldEnergy));
      buffer->addExtraInfo("opt-params", ExtraInfo(x));
      buffer->addExtraInfo("opt-ansatz", ExtraInfo(ansatzOps));
      
      // Persist XASM code to buffer
      auto xasm = xacc::getCompiler("xasm");
      auto xasmCode = xasm->translate(ansatz->operator()(x)); 
      buffer->addExtraInfo("opt-circuit", ExtraInfo(xasmCode));

    } else {
      xacc::info("ADAPT-" + subAlgo + " did not converge in " +
                 std::to_string(_maxIter) + " iterations.");
      return;
    }
  }
}

} // namespace algorithm
} // namespace xacc