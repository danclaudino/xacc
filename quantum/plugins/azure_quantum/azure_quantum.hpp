/*******************************************************************************
 * Copyright (c) 2023 UT-Battelle, LLC.
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
#ifndef XACC_AZURE_QUANTUM_ACCELERATOR_HPP_
#define XACC_AZURE_QUANTUM_ACCELERATOR_HPP_

#include "InstructionIterator.hpp"
#include "Accelerator.hpp"
#include <bitset>
#include <type_traits>
#include "azure_utils.hpp"

namespace xacc {
namespace quantum {

class AzureQuantumAccelerator : public Accelerator {
public:
  void initialize(const HeterogeneousMap &params = {}) override;
  void updateConfiguration(const HeterogeneousMap &config) override {
    if (config.keyExists<int>("shots")) {
      shots = config.get<int>("shots");
    }
    if (config.stringExists("backend")) {
      backend = config.getString("backend");
    }
    if (config.keyExists<bool>("http-verbose")) {
      //restClient->setVerbose(config.get<bool>("http-verbose"));
    }
  }

  const std::vector<std::string> configurationKeys() override {
    return {"shots", "backend"};
  }

  const std::string getSignature() override { return "azure-quantum"; }

  const std::string name() const override { return "azure-quantum"; }
  const std::string description() const override { return ""; }

  void execute(std::shared_ptr<AcceleratorBuffer> buffer,
               const std::shared_ptr<CompositeInstruction> circuit) override;

  void execute(std::shared_ptr<AcceleratorBuffer> buffer,
               const std::vector<std::shared_ptr<CompositeInstruction>>
                   circuits) override;

  bool isRemote() override { return true; }

  AzureQuantumAccelerator()
      : Accelerator() {}

  virtual ~AzureQuantumAccelerator() {}

private:
  void searchAPIKey(std::string &key);
  void findApiKeyInFile(std::string &key, const std::string &p);
  //std::map<std::string, std::string> generateRequestHeader(const std::string user_agent = "azure-quantum-cpp-demo-v0.1") const;

  std::string accessToken;
  std::string baseUrl;

  std::string job_name;
  std::string retrieve_job_id;

  int shots = 1024;
  std::string backend = "";
  bool initialized = false;

  // List of backend names:
  std::vector<std::string> available_backends;
  /*
  std::string post(const std::string &_url, const std::string &path,
                   const std::string &postStr,
                   std::map<std::string, std::string> headers = {});

  std::string get(const std::string &_url, const std::string &path,
                  std::map<std::string, std::string> headers =
                      std::map<std::string, std::string>{},
                  std::map<std::string, std::string> extraParams = {});
  */
};

} // namespace quantum
} // namespace xacc

#endif
