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
#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "xacc.hpp"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include "ObservableGrouping.hpp"

TEST(QWCTester, checkGrouping) {

  auto str = std::string(
      "Z0 + Z1 + Z2 + Z3 + Z0 Z1 + Z0 Z2 + Z0 Z3 + Z1 Z2 + Z1 Z3 + Z2 Z3 + "
      "X0 X1 Y2 Y3 + X0 Y1 Y2 X3 + Y0 X1 X2 Y3 + Y0 Y1 X2 X3");

  auto H = xacc::quantum::getObservable("pauli", str);
  auto qwc = xacc::getService<xacc::ObservableGrouping>("qwc");
  auto com = qwc->group(H);

  std::cout << com[0]->toString() << "\n";

  EXPECT_EQ("(1,0) X0 Y1 Y2 X3", com[0]->toString());
  EXPECT_EQ("(1,0) Y0 Y1 X2 X3", com[1]->toString());
  EXPECT_EQ("(1,0) Y0 X1 X2 Y3", com[2]->toString());
  EXPECT_EQ("(1,0) X0 X1 Y2 Y3", com[3]->toString());
  EXPECT_EQ("(1,0) Z1 Z2 + (1,0) Z0 Z3 + (1,0) Z1 Z3 + (1,0) Z2 Z3 + (1,0) Z0 "
            "Z1 + (1,0) Z0 Z2 + (1,0) Z3 + (1,0) Z2 + (1,0) Z0 + (1,0) Z1",
            com[4]->toString());
}

TEST(QWCTester, checkObserveGroups) {

  std::vector<std::shared_ptr<xacc::Observable>> ops;
  ops.push_back(xacc::quantum::getObservable(
      "pauli", std::string("Y1 X2 + Z0 Z1 Z2 + X0 X2")));

  auto qwc = xacc::getService<xacc::ObservableGrouping>("qwc");
  auto circuit = qwc->observeGroups(ops);

  std::string expected = R"expected(H q0
Swap q0,q1
CNOT q0,q2
S q0
CZ q1,q2
H q0
Measure q0
H q1
Measure q1
H q2
Measure q2
)expected";

  EXPECT_EQ(expected, circuit[0]->toString());
}

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}