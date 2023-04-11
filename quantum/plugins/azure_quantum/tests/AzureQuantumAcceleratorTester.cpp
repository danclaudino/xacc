#include "xacc.hpp"
#include <gtest/gtest.h>

TEST(AzureQuantumAcceleratorTester, checkInitialize) {
  xacc::set_verbose(true);
  auto accelerator = xacc::getAccelerator("azure-quantum");
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
}




int main(int argc, char **argv) {
  xacc::Initialize();
  //   xacc::set_verbose(true);
  ::testing::InitGoogleTest(&argc, argv);
  const auto result = RUN_ALL_TESTS();

  xacc::Finalize();

  return result;
}