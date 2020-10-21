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
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#include <gtest/gtest.h>
#include <utility>

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Optimizer.hpp"
#include "Observable.hpp"
#include "Algorithm.hpp"
#include "PauliOperator.hpp"
#include "xacc_observable.hpp"
#include "AlgorithmGradientStrategy.hpp"

using namespace xacc;
using namespace xacc::quantum;



TEST(RBMChemistryTester, checkH) {

  xacc::set_verbose(true);
  xacc::external::load_external_language_plugins();

  auto str = std::string("(-0.165606823582,-0)  2^ 1^ 2 1 + (0.120200490713,0)  2^ 0^ 0 2 + "
                          "(-0.0454063328691,-0)  0^ 3^ 2 1 + (0.168335986252,0)  1^ 0^ 0 1 + "
                          "(0.0454063328691,0)  2^ 1^ 3 0 + (0.168335986252,0)  0^ 1^ 1 0 + "
                          "(0.165606823582,0)  0^ 3^ 3 0 + (-0.0454063328691,-0)  3^ 0^ 1 2 + "
                          "(-0.0454063328691,-0)  2^ 3^ 0 1 + (-0.0454063328691,-0)  3^ 2^ 1 0 + "
                          "(0.165606823582,0)  2^ 1^ 1 2 + (-0.165606823582,-0)  0^ 3^ 0 3 + "
                          "(-0.479677813134,-0)  3^ 3 + (-0.0454063328691,-0)  2^ 1^ 0 3 + "
                          "(-0.174072892497,-0)  2^ 3^ 2 3 + (-0.0454063328691,-0)  0^ 1^ 2 3 + "
                          "(0.120200490713,0)  0^ 2^ 2 0 + (0.0454063328691,0)  0^ 1^ 3 2 + "
                          "(0.174072892497,0)  2^ 3^ 3 2 + (0.165606823582,0)  1^ 2^ 2 1 + "
                          "(-0.0454063328691,-0)  1^ 2^ 3 0 + (-0.120200490713,-0)  1^ 3^ 1 3 + "
                          "(0.120200490713,0)  1^ 3^ 3 1 + (-0.168335986252,-0)  0^ 1^ 0 1 + "
                          "(0.120200490713,0)  3^ 1^ 1 3 + (-0.120200490713,-0)  3^ 1^ 3 1 + "
                          "(0.0454063328691,0)  2^ 3^ 1 0 + (-1.2488468038,-0)  0^ 0 + "
                          "(0.0454063328691,0)  3^ 2^ 0 1 + (-0.168335986252,-0)  1^ 0^ 1 0 + "
                          "(0.165606823582,0)  3^ 0^ 0 3 + (-0.0454063328691,-0)  1^ 0^ 3 2 + "
                          "(0.0454063328691,0)  1^ 0^ 2 3 + (-1.2488468038,-0)  1^ 1 + "
                          "(0.0454063328691,0)  1^ 2^ 0 3 + (0.174072892497,0)  3^ 2^ 2 3 + "
                          "(-0.479677813134,-0)  2^ 2 + (-0.174072892497,-0)  3^ 2^ 3 2 + "
                          "(0.0454063328691,0)  3^ 0^ 2 1 + (-0.165606823582,-0)  3^ 0^ 3 0 + "
                          "(0.0454063328691,0)  0^ 3^ 1 2 + (-0.165606823582,-0)  1^ 2^ 1 2 + "
                          "(-0.120200490713,-0)  0^ 2^ 0 2 + (-0.120200490713,-0)  2^ 0^ 2 0 + (0.7080240981,0)");


  auto H = xacc::quantum::getObservable("fermion", str);
  auto rbm = xacc::getService<xacc::Algorithm>("rbm-chemistry");
  //auto acc = xacc::getAccelerator("dwave-ocean-sdk", {{"backend", "DW_2000Q_6"}, {"shots", 1000}});
   // auto acc = xacc::getAccelerator("dwave:DW_2000Q_6");
    auto acc = xacc::getAccelerator("dwave-neal", {{"shots", 1000}});
  EXPECT_TRUE(rbm->initialize({std::make_pair("hidden-layers", 4), std::make_pair("accelerator", acc),
                                std::make_pair("observable", H)}));

  auto buffer = xacc::qalloc(4);
  rbm->execute(buffer);

  xacc::external::unload_external_language_plugins();

}

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  auto ret = RUN_ALL_TESTS();
  xacc::Finalize();
  return ret;
}