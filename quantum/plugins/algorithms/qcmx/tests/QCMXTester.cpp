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
 ******************************************************************************/
#include <gtest/gtest.h>
#include <memory>
#include <utility>
#include <vector>
#include <iomanip>

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"

using namespace xacc;
using namespace xacc::quantum;

TEST(QCMXTester, checkSimple) {

  auto H = xacc::quantum::getObservable(
      "pauli", std::string("0.2976 + 0.3593 Z0 - 0.4826 Z1 + 0.5818 Z0 Z1 + "
                           "0.0896 X0 X1 + 0.0896 Y0 Y1"));

  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  ansatz->addInstruction(provider->createInstruction("X", {(size_t)0}));

  auto acc = xacc::getAccelerator("qpp");

  auto qcmx = xacc::getService<xacc::Algorithm>("qcmx");
  qcmx->initialize({{"accelerator", acc},
                    {"observable", H},
                    {"ansatz", ansatz},
                    {"cmx-order", 2}});

  auto buffer = qalloc(2);
  qcmx->execute(buffer);
  auto map =
      buffer->getInformation("energies").as<std::map<std::string, double>>();

std::cout <<  map["PDS(2)"] << "\n";
  EXPECT_NEAR(-1.14499, map["PDS(2)"], 1e-3);
  EXPECT_NEAR(-1.14517, map["CMX(2)"], 1e-3);
  EXPECT_NEAR(-1.14517, map["Knowles(2)"], 1e-3);
  EXPECT_NEAR(-1.14499, map["Soldatov"], 1e-3);
}

/*
TEST(QCMXTester, checkHA) {

  auto H = xacc::quantum::getObservable(
      "pauli", std::string("1.5 + 0.5 Z1 - Z0 Z1"));

  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  ansatz->addVariables({"x0", "x1"});
  ansatz->addInstruction(provider->createInstruction("Rx", {0}, {"x0"}));
  ansatz->addInstruction(provider->createInstruction("CNOT", {0, 1}));
  ansatz->addInstruction(provider->createInstruction("Ry", {1}, {"-x1"}));
  ansatz->addInstruction(provider->createInstruction("CNOT", {0, 1}));
  ansatz->addInstruction(provider->createInstruction("Ry", {1}, {"x1"}));

  //std::cout << ansatz->toString() << "\n";

  auto acc = xacc::getAccelerator("qsim");

  auto qcmx = xacc::getService<xacc::Algorithm>("pds-vqs");

  auto optimizer = xacc::getOptimizer("nlopt", {{"step", 0.05}, {"maxeval", 10}, {"algorithm", "l-bfgs"}, {"initial-parameters", std::vector<double>{xacc::constants::pi / 2.0, xacc::constants::pi / 4.0} }});

  qcmx->initialize({{"accelerator", acc},
                    {"observable", H},
                    {"ansatz", ansatz},
                    {"metric", "ITE"},
                    {"optimizer", optimizer},
                    {"cmx-order", 2}});

  auto buffer = qalloc(2);
  qcmx->execute(buffer);
  auto map =
      buffer->getInformation("opt-val").as<double>();

    std::cout <<  map << "\n";
  //EXPECT_NEAR(-1.14499, map["PDS(2)"], 1e-3);

}
*/

TEST(QCMXTester, checkPDS_VQS) {

  xacc::set_verbose(true);

  auto H = xacc::quantum::getObservable(
      "pauli", std::string("0.4 Z0 + 0.4 Z1 + 0.2 X0 X1"));

  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  ansatz->addVariables({"x0", "x1", "x2", "x3"});
  ansatz->addInstruction(provider->createInstruction("Ry", {0}, {"x0"}));
  ansatz->addInstruction(provider->createInstruction("Ry", {1}, {"x1"}));
  ansatz->addInstruction(provider->createInstruction("CNOT", {0, 1}));
  ansatz->addInstruction(provider->createInstruction("Ry", {1}, {"x2"}));
  ansatz->addInstruction(provider->createInstruction("Ry", {0}, {"x3"}));

  auto acc = xacc::getAccelerator("qsim");
  auto optimizer = xacc::getOptimizer("nlopt", {{"maxeval", 10}, {"algorithm", "l-bfgs"}, {"initial-parameters", std::vector<double>{7 * xacc::constants::pi / 16.0, xacc::constants::pi - 0.0001, 0.0 , 0.0}}});


  auto qcmx = xacc::getService<xacc::Algorithm>("pds-vqs");
  qcmx->initialize({{"accelerator", acc},
                    {"observable", H},
                    {"ansatz", ansatz},
                    {"optimizer", optimizer},
                    {"cmx-order", 2}});

  auto buffer = qalloc(2);
  qcmx->execute(buffer);

}


TEST(QCMXTester, checkADAPT_PDS_VQS) {

  xacc::set_verbose(true);
  //xacc::logToFile(true, "./log");
  auto acc = xacc::getAccelerator("qsim");

  auto optimizer = xacc::getOptimizer("nlopt", {std::make_pair("nlopt-optimizer", "l-bfgs")});
  auto adapt_pds_vqs = xacc::getService<xacc::Algorithm>("pds-vqs");
  int nElectrons = 2;
  auto pool = "singlet-adapted-uccsd";

  auto str = std::string("(-0.165606823582,-0)  1^ 2^ 1 2 + (0.120200490713,0)  1^ 0^ 0 1 + "
                          "(-0.0454063328691,-0)  0^ 3^ 1 2 + (0.168335986252,0)  2^ 0^ 0 2 + "
                          "(0.0454063328691,0)  1^ 2^ 3 0 + (0.168335986252,0)  0^ 2^ 2 0 + "
                          "(0.165606823582,0)  0^ 3^ 3 0 + (-0.0454063328691,-0)  3^ 0^ 2 1 + "
                          "(-0.0454063328691,-0)  1^ 3^ 0 2 + (-0.0454063328691,-0)  3^ 1^ 2 0 + "
                          "(0.165606823582,0)  1^ 2^ 2 1 + (-0.165606823582,-0)  0^ 3^ 0 3 + "
                          "(-0.479677813134,-0)  3^ 3 + (-0.0454063328691,-0)  1^ 2^ 0 3 + "
                          "(-0.174072892497,-0)  1^ 3^ 1 3 + (-0.0454063328691,-0)  0^ 2^ 1 3 + "
                          "(0.120200490713,0)  0^ 1^ 1 0 + (0.0454063328691,0)  0^ 2^ 3 1 + "
                          "(0.174072892497,0)  1^ 3^ 3 1 + (0.165606823582,0)  2^ 1^ 1 2 + "
                          "(-0.0454063328691,-0)  2^ 1^ 3 0 + (-0.120200490713,-0)  2^ 3^ 2 3 + "
                          "(0.120200490713,0)  2^ 3^ 3 2 + (-0.168335986252,-0)  0^ 2^ 0 2 + "
                          "(0.120200490713,0)  3^ 2^ 2 3 + (-0.120200490713,-0)  3^ 2^ 3 2 + "
                          "(0.0454063328691,0)  1^ 3^ 2 0 + (-1.2488468038,-0)  0^ 0 + "
                          "(0.0454063328691,0)  3^ 1^ 0 2 + (-0.168335986252,-0)  2^ 0^ 2 0 + "
                          "(0.165606823582,0)  3^ 0^ 0 3 + (-0.0454063328691,-0)  2^ 0^ 3 1 + "
                          "(0.0454063328691,0)  2^ 0^ 1 3 + (-1.2488468038,-0)  2^ 2 + "
                          "(0.0454063328691,0)  2^ 1^ 0 3 + (0.174072892497,0)  3^ 1^ 1 3 + "
                          "(-0.479677813134,-0)  1^ 1 + (-0.174072892497,-0)  3^ 1^ 3 1 + "
                          "(0.0454063328691,0)  3^ 0^ 1 2 + (-0.165606823582,-0)  3^ 0^ 3 0 + "
                          "(0.0454063328691,0)  0^ 3^ 2 1 + (-0.165606823582,-0)  2^ 1^ 2 1 + "
                          "(-0.120200490713,-0)  0^ 1^ 0 1 + (-0.120200490713,-0)  1^ 0^ 1 0 + (0.7080240981,0)");

  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  ansatz->addInstruction(provider->createInstruction("X", {0}));
  ansatz->addInstruction(provider->createInstruction("X", {2}));

  auto H = xacc::quantum::getObservable("fermion", str);

  EXPECT_TRUE(adapt_pds_vqs->initialize({std::make_pair("accelerator",acc),
                                std::make_pair("observable", H),
                                std::make_pair("optimizer", optimizer),
                                std::make_pair("pool", pool),
                                std::make_pair("ansatz", ansatz),
                                std::make_pair("adapt", true),
                                std::make_pair("cmx-order", 2),
                                std::make_pair("n-electrons", nElectrons)
                                }));

  auto buffer = xacc::qalloc(4);
  adapt_pds_vqs->execute(buffer);
  //EXPECT_NEAR(-1.13717, buffer_vqe->getInformation("opt-val").as<double>(), 1e-4);
}

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}