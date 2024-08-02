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
 *   Thien Nguyen - initial API and implementation
 *******************************************************************************/

#include "Accelerator.hpp"
#include "xacc.hpp"
#include <gtest/gtest.h>
#include "xacc_observable.hpp"

TEST(NewAerTester, checkSimple) {
  auto accelerator =
      xacc::getAccelerator("aer", {std::make_pair("shots", 2048)});
  auto xasmCompiler = xacc::getCompiler("xasm");
  auto ir = xasmCompiler->compile(R"(__qpu__ void bell(qbit q) {
      H(q[0]);
      CX(q[0], q[1]);
      Measure(q[0]);
      Measure(q[1]);
    })",
                                  accelerator);

  auto program = ir->getComposite("bell");

  auto buffer = xacc::qalloc(2);
  accelerator->execute(buffer, program);

  buffer->print(std::cout);

  accelerator->updateConfiguration({std::make_pair("sim-type", "statevector")});

  buffer->resetBuffer();
  accelerator->execute(buffer, program);
  buffer->print(std::cout);
}


TEST(NewAerTester, checkDeuteron) {
  auto accelerator =
      xacc::getAccelerator("aer", {std::make_pair("sim-type", "statevector")});
  auto xasmCompiler = xacc::getCompiler("xasm");
  auto ir = xasmCompiler->compile(R"(__qpu__ void ansatz(qbit q, double t) {
      X(q[0]);
      Ry(q[1], t);
      CX(q[1], q[0]);
      H(q[0]);
      H(q[1]);
      Measure(q[0]);
      Measure(q[1]);
    })",
                                  accelerator);

  auto program = ir->getComposite("ansatz");
  // Expected results from deuteron_2qbit_xasm_X0X1
  const std::vector<double> expectedResults{
      0.0,       -0.324699, -0.614213, -0.837166, -0.9694,
      -0.996584, -0.915773, -0.735724, -0.475947, -0.164595,
      0.164595,  0.475947,  0.735724,  0.915773,  0.996584,
      0.9694,    0.837166,  0.614213,  0.324699,  0.0};

  const auto angles =
      xacc::linspace(-xacc::constants::pi, xacc::constants::pi, 20);
  for (size_t i = 0; i < angles.size(); ++i) {
    auto buffer = xacc::qalloc(2);
    auto evaled = program->operator()({angles[i]});
    accelerator->execute(buffer, evaled);
    EXPECT_NEAR(buffer->getExpectationValueZ(), expectedResults[i], 1e-6);
  }
}

#ifndef QIREE_BUILD
TEST(NewAerTester, checkDeuteronVqeH2) {
  auto accelerator =
      xacc::getAccelerator("aer", {std::make_pair("sim-type", "statevector")});

  // Create the N=2 deuteron Hamiltonian
  auto H_N_2 = xacc::quantum::getObservable(
      "pauli", std::string("5.907 - 2.1433 X0X1 "
                           "- 2.1433 Y0Y1"
                           "+ .21829 Z0 - 6.125 Z1"));

  auto optimizer = xacc::getOptimizer("nlopt");
  xacc::qasm(R"(
        .compiler xasm
        .circuit deuteron_ansatz
        .parameters theta
        .qbit q
        X(q[0]);
        Ry(q[1], theta);
        CNOT(q[1],q[0]);
    )");
  auto ansatz = xacc::getCompiled("deuteron_ansatz");

  // Get the VQE Algorithm and initialize it
  auto vqe = xacc::getAlgorithm("vqe");
  vqe->initialize({std::make_pair("ansatz", ansatz),
                   std::make_pair("observable", H_N_2),
                   std::make_pair("accelerator", accelerator),
                   std::make_pair("optimizer", optimizer)});

  // Allocate some qubits and execute
  auto buffer = xacc::qalloc(2);
  vqe->execute(buffer);

  // Expected result: -1.74886
  EXPECT_NEAR((*buffer)["opt-val"].as<double>(), -1.74886, 1e-4);
}

TEST(NewAerTester, testDeuteronVqeH3) {
  auto accelerator = xacc::getAccelerator("aer", {{"sim-type", "statevector"}});

  // Create the N=3 deuteron Hamiltonian
  auto H_N_3 = xacc::quantum::getObservable(
      "pauli",
      std::string("5.907 - 2.1433 X0X1 - 2.1433 Y0Y1 + .21829 Z0 - 6.125 Z1 + "
                  "9.625 - 9.625 Z2 - 3.91 X1 X2 - 3.91 Y1 Y2"));

  auto optimizer = xacc::getOptimizer("nlopt");

  // JIT map Quil QASM Ansatz to IR
  xacc::qasm(R"(
        .compiler xasm
        .circuit deuteron_ansatz_h3
        .parameters t0, t1
        .qbit q
        X(q[0]);
        exp_i_theta(q, t1, {{"pauli", "X0 Y1 - Y0 X1"}});
        exp_i_theta(q, t0, {{"pauli", "X0 Z1 Y2 - X2 Z1 Y0"}});
    )");
  auto ansatz = xacc::getCompiled("deuteron_ansatz_h3");

  // Get the VQE Algorithm and initialize it
  auto vqe = xacc::getAlgorithm("vqe");
  vqe->initialize({std::make_pair("ansatz", ansatz),
                   std::make_pair("observable", H_N_3),
                   std::make_pair("accelerator", accelerator),
                   std::make_pair("optimizer", optimizer)});

  // Allocate some qubits and execute
  auto buffer = xacc::qalloc(3);
  vqe->execute(buffer);

  // Expected result: -2.04482
  EXPECT_NEAR((*buffer)["opt-val"].as<double>(), -2.04482, 1e-4);
}

#endif

TEST(NewAerTester, testExecutionInfoStateVec) {
  auto accelerator =
      xacc::getAccelerator("aer", {std::make_pair("sim-type", "statevector")});

  xacc::qasm(R"(
        .compiler xasm
        .circuit test_bell_exe
        .qbit q
        H(q[0]);
        CNOT(q[0],q[1]);
    )");
  auto bell = xacc::getCompiled("test_bell_exe");

  // Allocate some qubits and execute
  auto buffer = xacc::qalloc(2);
  accelerator->execute(buffer, bell);

  auto exeInfo = accelerator->getExecutionInfo();
  EXPECT_GT(exeInfo.size(), 0);
  auto waveFn =
      accelerator->getExecutionInfo<xacc::ExecutionInfo::WaveFuncPtrType>(
          xacc::ExecutionInfo::WaveFuncKey);
  for (const auto &elem : *waveFn) {
    std::cout << elem << "\n";
  }
  // 2 qubits => 4 elements
  EXPECT_EQ(waveFn->size(), 4);
  EXPECT_NEAR(std::abs((*waveFn)[0] - 1.0 / std::sqrt(2.0)), 0.0, 1e-9);
  EXPECT_NEAR(std::abs((*waveFn)[3] - 1.0 / std::sqrt(2.0)), 0.0, 1e-9);
}

TEST(NewAerTester, checkInitialState1) {
  const std::vector<std::complex<double>> initial_states{0.0, 1.0};
  auto accelerator = xacc::getAccelerator(
      "aer", {{"sim-type", "statevector"}, {"initial_state", initial_states}});
  auto xasmCompiler = xacc::getCompiler("xasm");
  auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
      H(q[0]);
    })",
                                  accelerator);

  auto program = ir->getComposite("test");

  auto buffer = xacc::qalloc(1);
  accelerator->execute(buffer, program);

  auto exeInfo = accelerator->getExecutionInfo();
  EXPECT_GT(exeInfo.size(), 0);
  auto waveFn =
      accelerator->getExecutionInfo<xacc::ExecutionInfo::WaveFuncPtrType>(
          xacc::ExecutionInfo::WaveFuncKey);
  for (const auto &elem : *waveFn) {
    std::cout << elem << "\n";
  }
  EXPECT_EQ(waveFn->size(), 2);
  // Expect |-> state: |0> - |1> since we set the initial state to |1>
  EXPECT_NEAR(std::abs((*waveFn)[0] - 1.0 / std::sqrt(2.0)), 0.0, 1e-9);
  EXPECT_NEAR(std::abs((*waveFn)[1] + 1.0 / std::sqrt(2.0)), 0.0, 1e-9);
}

int main(int argc, char **argv) {
  xacc::Initialize();

  ::testing::InitGoogleTest(&argc, argv);
  const auto result = RUN_ALL_TESTS();

  xacc::Finalize();

  return result;
}
