#include "xacc.hpp"
#include "xacc_observable.hpp"

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);

  // Get the desired Accelerator and Optimizer
  auto qpu = xacc::getAccelerator("tnqvm");
  auto optimizer = xacc::getOptimizer("nlopt");

  // Create the N=3 deuteron Hamiltonian
  auto H2 = xacc::quantum::getObservable(
     "pauli",
     std::string("(0.174073,0) Z2 Z3 + (0.1202,0) Z1 Z3 + (0.165607,0) Z1 Z2 + "
        "(0.165607,0) Z0 Z3 + (0.1202,0) Z0 Z2 + (-0.0454063,0) Y0 Y1 X2 X3 + "
        "(-0.220041,0) Z3 + (-0.106477,0) + (0.17028,0) Z0 + (-0.220041,0) Z2 "
        "+ (0.17028,0) Z1 + (-0.0454063,0) X0 X1 Y2 Y3 + (0.0454063,0) X0 Y1 "
        "Y2 X3 + (0.168336,0) Z0 Z1 + (0.0454063,0) Y0 X1 X2 Y3"));

  // Get the VQE Algorithm and initialize it
  auto adapt_vqe = xacc::getAlgorithm("adapt-vqe");
  adapt_vqe->initialize({
        std::make_pair("observable", H2),
        std::make_pair("accelerator", qpu),
        std::make_pair("optimizer", optimizer),
        std::make_pair("pool", "uccsd"),
        std::make_pair("nElectrons", 2)
        });

  // Allocate some qubits and execute
  auto buffer = xacc::qalloc(4);
  adapt_vqe->execute(buffer);

  // Get the result
  //auto energy =
        //(*buffer)["opt-val"].as<double>();
  //std::cout << "Energy: " << energy << "\n";
  xacc::Finalize();
  return 0;
}