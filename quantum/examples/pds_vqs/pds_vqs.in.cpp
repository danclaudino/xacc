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
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include <fstream>

int main(int argc, char **argv) {
  xacc::set_verbose(true);
  // Initialize and load the python plugins
  xacc::Initialize(argc, argv);

  // Process the input arguments
  std::vector<std::string> arguments(argv + 1, argv + argc);
  std::string geometry, opt = "cobyla";
  int order = 2;
  double tol = 1.0e-6;
  for (int i = 0; i < arguments.size(); i++) {
    if (arguments[i] == "--pool") {
      order = std::stoi(arguments[i + 1]);
    }
    if (arguments[i] == "--geometry") {
      geometry = arguments[i + 1];
    }
    if (arguments[i] == "--tol") {
      tol = std::stod(arguments[i + 1]);
    }
    if (arguments[i] == "--optimizer") {
      opt = arguments[i + 1];
    }
  }

  auto accelerator = xacc::getAccelerator("qsim");

  std::ifstream inFile;
  inFile.open("@CMAKE_SOURCE_DIR@/quantum/examples/pds_vqs/h4_" + geometry + ".txt");
  std::stringstream strStream;
  strStream << inFile.rdbuf();
  auto H = xacc::quantum::getObservable("fermion", strStream.str());

  auto q = xacc::qalloc(H->nBits());

  auto optimizer = xacc::getOptimizer("nlopt", {{"algorithm", opt}, {"ftol", tol}});

  auto tmp = xacc::getService<xacc::Instruction>("uccsd");
  auto uccsd = std::dynamic_pointer_cast<xacc::CompositeInstruction>(tmp);
	uccsd->expand({std::make_pair("ne", 4), std::make_pair("nq", 8)});
  auto provider = xacc::getService<xacc::IRProvider>("quantum");
  auto ansatz = provider->createComposite("initial-state");
  for (std::size_t i : {0, 1, 4, 5}) {
    ansatz->addInstruction(provider->createInstruction("X", {i}));
  }
  ansatz->addVariables(uccsd->getVariables());
  for (auto & inst : uccsd->getInstructions()) ansatz->addInstruction(inst);

  auto pds_vqs = xacc::getAlgorithm("pds-vqs", {{"observable", H},
                                                {"accelerator", accelerator},
                                                {"ansatz", uccsd},
                                                {"adapt", false},
                                                {"cmx-order", order},
                                                {"optimizer", optimizer}});

  pds_vqs->execute(q);
  std::cout << "Optimum PDS(" << order << ")-VQS energy: " << q->getInformation("opt-val").as<double>() << "\n\n";
  std::cout << "Optimum parameters: ";
  for (auto x : q->getInformation("opt-params").as<std::vector<double>>()) {
    std::cout << x << "  ";
  }
  //
  // Finalize
  xacc::Finalize();
}
