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
 *   Daniel Claudino - Update to Qiskit Runtime (Qiskit 1.x)
 *******************************************************************************/
#include "IBMAccelerator.hpp"
#include <cctype>
#include <fstream>
#include <iostream>
#include <string>

#include "Properties.hpp"
#include "QObjectExperimentVisitor.hpp"
#include "CountGatesOfTypeVisitor.hpp"
#include "Utils.hpp"
#include "cpr/response.h"

#ifndef REMOTE_DISABLED
#include <cpr/cpr.h>
#endif

#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Scheduler.hpp"
#include "QObject.hpp"


#include <regex>
#include <thread>
#include <cassert>
namespace xacc {
namespace quantum {
const std::string IBMAccelerator::IBM_AUTH_URL =
    "https://auth.quantum-computing.ibm.com";
const std::string IBMAccelerator::IBM_API_URL =
"https://api.quantum-computing.ibm.com";

const std::string IBMAccelerator::IBM_TRANSPILER_URL = "https://cloud-transpiler.quantum.ibm.com";


const std::string IBMAccelerator::DEFAULT_IBM_BACKEND = "ibmq_qasm_simulator";
const std::string IBMAccelerator::IBM_LOGIN_PATH = "/api/users/loginWithToken";

std::string hex_string_to_binary_string(std::string hex) {
  return integral_to_binary_string((int)strtol(hex.c_str(), NULL, 0));
}

bool caseInsensitiveCompare(const std::string& str1, const std::string& str2) {
    return str1.size() == str2.size() &&
           std::equal(str1.begin(), str1.end(), str2.begin(), [](unsigned char c1, unsigned char c2) {
               return std::tolower(c1) == std::tolower(c2);
           });
}

bool hasMidCircuitMeasurement(
    const std::shared_ptr<CompositeInstruction> &in_circuit) {
  InstructionIterator it(in_circuit);
  bool measureEncountered = false;
  while (it.hasNext()) {
    auto nextInst = it.next();
    if (nextInst->isEnabled()) {
      if (nextInst->name() == "Measure") {
        // Flag that we have seen a Measure gate.
        measureEncountered = true;
      }

      // We have seen a Measure gate but this one is not another Measure gate.
      if (measureEncountered && nextInst->name() != "Measure") {
        // This circuit has mid-circuit measurement
        return true;
      }
    }
  }
  return false;
}

bool IBMAccelerator::verifyJobsLimit(std::string& curr_backend) {

  // Get backend jobs limit
  std::string getJobsLimitPath = "/api/Network/" + hub + "/Groups/" + group +
                        "/Projects/" + project + "/devices/" + curr_backend +
                        "/jobsLimit";
  auto backend_jobslimit_response = get(IBM_API_URL, getJobsLimitPath, {},
        {std::make_pair("version", "1"),
        std::make_pair("access_token", currentApiToken)});
  auto backend_jobslimit_response_json = json::parse(backend_jobslimit_response);

  auto maximumJobs = backend_jobslimit_response_json["maximumJobs"].get<int>();
  auto runningJobs = backend_jobslimit_response_json["runningJobs"].get<int>();

  return maximumJobs < 0 || runningJobs < maximumJobs;
}

void IBMAccelerator::processBackendCandidate(const nlohmann::json& backend_json) {

  // First of all filter by count of qubits
  if (requested_n_qubits > 0) {
    if( backend_json.count("n_qubits") ) {
      int nqubits = backend_json["n_qubits"].get<int>();
      if(nqubits < requested_n_qubits) {
        return;
      }
    } else {
      return;
    }
  }
  std::string curr_backend = backend_json["backend_name"].get<std::string>();

// Get current backend status
  const std::string path("/runtime/backends/" + curr_backend + "/status");
  auto status_json = json::parse(get(IBM_API_URL, path, headers, {{"name", curr_backend}}));

  auto queue_length = status_json["length_queue"].get<int>();
  auto state = status_json["state"].get<bool>();

  if (state && (backendQueueLength < 0 || backendQueueLength > queue_length)) {
    if (filterByJobsLimit && !verifyJobsLimit(curr_backend)) {
      return;
    }
    backendQueueLength = queue_length;
    auto old_backend = backend;
    backend = curr_backend;
    //availableBackends.clear();
    availableBackends.push_back(curr_backend);
  }
}

void IBMAccelerator::selectBackend(std::vector<std::string>& all_available_backends) {
  bool lowest_queue_backend = false;
  if (backend == "lowest-queue-count") {
    lowest_queue_backend = true;
  }

  for (std::string b : backends_root["devices"]) {

      const std::string path("/runtime/backends/" + b + "/configuration");
      auto backend_json = json::parse(get(IBM_API_URL, path, headers, {{"name", b}}));
  
    // Simple case: select by backend_name
    if (!lowest_queue_backend) {

      if (b == backend) {
        //const std::string path("/runtime/backends/" + backend + "/status");
        //auto backend_json = json::parse(get(IBM_API_URL, path, headers, {{"name", b}}));

        if(caseInsensitiveCompare(backend_json["status"].get<std::string>(), "active")) {
          //availableBackends.insert(std::make_pair(backend, backend_json));
          availableBackends.push_back(b);
        } else {
          xacc::error(backend + "not active at the moment.");
        }
      }
    } else {
      // Select backend by job queue size and by parameters (optional)
      //const std::string path("/runtime/backends/" + b + "/configuration");
      //auto backend_json = json::parse(get(IBM_API_URL, path, headers, {{"name", b}})); 
      processBackendCandidate(backend_json);
    }
    //all_available_backends.push_back(b["backend_name"].get<std::string>());
    all_available_backends.push_back(b);
    
  }
  if (lowest_queue_backend) {
    xacc::info("Backend with lowest queue count: " + backend);
  }
}


void IBMAccelerator::initialize(const HeterogeneousMap &params) {
  if (!initialized) {
    std::string apiKey = "";
    std::string getBackendPath = "/api/Backends?access_token=";
    std::string tokenParam = "{\"apiToken\": \"";
    std::string getBackendPropertiesPath;

    // Set backend, shots, etc.
    // and get the apikey, hub, group, and project
    updateConfiguration(params);
    searchAPIKey(apiKey, hub, group, project);
    currentApiToken = apiKey;

    // We should have backend set by now

    if (hub.empty() && group.empty() && project.empty()) {
      // Fallback to public API credentials
      hub = "ibm-q";
      group = "open";
      project = "main";
    }

    headers = {
        {"Accept", "application/json"},
        {"Content-Type", "application/json"},
        {"Authorization", "Bearer " + apiKey}
        };

    auto getUser = get(IBM_API_URL, "/runtime/users/me", headers);
    auto provider = json::parse(getUser)["instances"][0]["name"].get<std::string>();

    backends_root = json::parse(get(IBM_API_URL, "/runtime/backends", headers, {{"provider", provider}}));
    
    std::vector<std::string> your_available_backends;
    selectBackend(your_available_backends);

    getBackendPropertiesPath = "/runtime/backends/" + backend + "/properties";
    // Get current backend properties
    auto backend_props_response =
        get(IBM_API_URL, getBackendPropertiesPath, headers, {{"name", backend}});

    xacc::info("Backend property:\n" + backend_props_response);
    auto props = json::parse(backend_props_response);
    backendProperties.insert({backend, props});

    if (!xacc::container::contains(your_available_backends, backend)) {
      std::stringstream error_ss;
      error_ss << "IBM Initialization Error:\n";
      error_ss << "Hub: " << hub << "\n";
      error_ss << "Group: " << group << "\n";
      error_ss << "Project: " << project << "\n";
      error_ss << "The requested backend (" << backend
               << ") is not available in this allocation.";
      error_ss << "\n\nAvailable backends are:\n";
      for (int i = 0; i < your_available_backends.size(); i++) {
        error_ss << your_available_backends[i] << ( i < your_available_backends.size()-1 ? ", " : "");
        if (i % 4 == 0) error_ss << "\n";
      }
      xacc::error(error_ss.str());
    }


    const std::string path("/runtime/backends/" + backend + "/configuration");
    chosenBackend = json::parse(get(IBM_API_URL, path, headers, {{"name", backend}}));
    xacc::info("Backend config:\n" + chosenBackend.dump());

    multi_meas_enabled = chosenBackend.value("multi_meas_enabled", false);

    auto defaultsPath = "/runtime/backends/" + backend + "/defaults";
    defaults_response = get(IBM_API_URL, defaultsPath, headers);
    xacc::info("Backend default:\n" + defaults_response);

    // check if cloud transpiler is working in case it will be used
    if (useCloudTranspiler) {
      auto checkTranspilerHealth = get(IBM_TRANSPILER_URL, "/health", headers);
      auto transpilerHealth = json::parse(checkTranspilerHealth)["message"].get<std::string>();
      if (transpilerHealth != "Health OK") {
        xacc::warning("IBM cloud transpiler is not working. Unless your program is already transpiled, it will not run.");
      }
    }

    initialized = true;
  }

}

void IBMAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> circuit) {
  execute(buffer, std::vector<std::shared_ptr<CompositeInstruction>>{circuit});
}

std::string QasmQObjGenerator::getQObjJsonStr(
    std::vector<std::shared_ptr<CompositeInstruction>> circuits,
    const int &shots, const nlohmann::json &backend,
    const std::string getBackendPropsResponse,
    std::vector<std::pair<int, int>> &connectivity,
    const nlohmann::json &backendDefaults) {

    std::cout << "QasmQObjGenerato\n";
  // Create a QObj
  xacc::ibm::QObject qobj;
  qobj.set_qobj_id("xacc-qobj-id");
  qobj.set_schema_version("1.1.0");
  qobj.set_type("QASM");
  qobj.set_header(QObjectHeader());

  const auto basis_gates =
      backend["basis_gates"].get<std::vector<std::string>>();

      std::cout << basis_gates << "\n";
  // If the gate set has "u3" -> old gateset.
  const auto gateSet = (xacc::container::contains(basis_gates, "cz"))
                           ? QObjectExperimentVisitor::GateSet::U_CX
                           : QObjectExperimentVisitor::GateSet::RZ_SX_CX;
  // Create the Experiments
  std::vector<xacc::ibm::Experiment> experiments;
  int maxMemSlots = 0;
  for (auto &kernel : circuits) {

    auto visitor = std::make_shared<QObjectExperimentVisitor>(
        kernel->name(), backend["n_qubits"].get<int>(), gateSet);

    InstructionIterator it(kernel);
    int memSlots = 0;
    while (it.hasNext()) {
      auto nextInst = it.next();
      if (nextInst->isEnabled()) {
        nextInst->accept(visitor);
      }
    }

    // After calling getExperiment, maxMemorySlots should be
    // maxClassicalBit + 1
    auto experiment = visitor->getExperiment();
    experiments.push_back(experiment);
    if (visitor->maxMemorySlots > maxMemSlots) {
      maxMemSlots = visitor->maxMemorySlots;
    }

    // Ensure that this Experiment maps onto the
    // hardware connectivity. Before now, we assume
    // any IRTransformations have been run.
    for (auto &inst : experiment.get_instructions()) {
      if (inst.get_name() == "cx") {

        std::set<int> connection;
        connection.insert(inst.get_qubits()[0]);
        connection.insert(inst.get_qubits()[1]);

        auto it =
            std::find_if(connectivity.begin(), connectivity.end(),
                         [&connection](const std::pair<int, int> &element) {
                           return std::set<int>{element.first, element.second} == connection;
                         });
        if (it == std::end(connectivity)) {
          std::stringstream ss;
          ss << "Invalid logical program connectivity, no connection between "
             << inst.get_qubits();
          xacc::error(ss.str());
        }
      }
    }
  }

  // Create the QObj Config
  xacc::ibm::QObjectConfig config;
  config.set_shots(shots);
  config.set_memory(false);
  config.set_meas_level(2);
  config.set_memory_slots(maxMemSlots);
  config.set_meas_return("avg");
  config.set_memory_slot_size(100);
  config.set_n_qubits(backend["n_qubits"].get<int>());

  // Add the experiments and config
  qobj.set_experiments(experiments);
  qobj.set_config(config);

  // Set the Backend
  xacc::ibm::Backend bkend;
  bkend.set_name(backend["backend_name"].get<std::string>());

  xacc::ibm::QObjectHeader qobjHeader;
  qobjHeader.set_backend_version("1.0.0");
  qobjHeader.set_backend_name(backend["backend_name"].get<std::string>());
  qobj.set_header(qobjHeader);

  // Create the JSON String to send
  nlohmann::json j = qobj;

  return j.dump();
}

void IBMAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>> circuits) {

    auto qasmCompiler = xacc::getCompiler("staq");
    std::vector<std::string> qasmCircuits;
    for (auto c : circuits) {
      qasmCircuits.push_back(qasmCompiler->translate(c));
    }

    // transpilation
     json transpilerBodyParams = {
        {"qasm_circuits", qasmCircuits}
    };

    json transpilerQueryParams = {{"backend", backend}};


    auto transpilerResponse = post(IBM_TRANSPILER_URL, "/transpile", transpilerBodyParams.dump(), headers, transpilerQueryParams.dump());

    auto task_id = json::parse(transpilerResponse)["task_id"].get<std::string>();

  int dots = 1;
  std::string state;
  while (true) {
    
    auto get_job_status = get(IBM_TRANSPILER_URL, "/transpile/" + task_id, headers);
    json get_job_status_json;

    try {
      get_job_status_json = json::parse(get_job_status);
    } catch (json::exception const &) {
      std::cout << "Failed to parse response '" << get_job_status
                << "' from job status " << task_id << std::endl;
      throw;
    }
    state = get_job_status_json["state"].get<std::string>();

    if (state.find("SUCCESS") == std::string::npos && state.find("PENDING") == std::string::npos) {
      xacc::error("Transpiler failed: " + get_job_status_json.dump(4));
    }

    if (state == "SUCCESS") {
      std::cout << std::endl;
      break;
    }

    if (dots > 4)
      dots = 1;
    std::stringstream ss;
    ss << "\033[0;32m"
        << "IBM Transpiler Job "
        << "\033[0;36m" << task_id << "\033[0;32m"
        << " Status: " << state;
    for (int i = 0; i < dots; i++)
      ss << '.';
    dots++;
    ss << "\033[0m";
    std::cout << '\r' << ss.str() << std::setw(20) << std::setfill(' ')
              << std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }

auto transpiled = json::parse(get(IBM_TRANSPILER_URL, "/transpile/" + task_id, headers));

    nlohmann::json body = {
        {"program_id", "sampler"},
        {"hub", hub},
        {"group", group},
        {"project", project},
        {"backend", backend},
        {"params", {
            {"pubs", json::array()},
            {"support_qiskit", false},
            {"version", 2}
        }}
    };

  for (const auto & c : transpiled["result"]) {
    body["params"]["pubs"].push_back(
      {c["qasm"], json::array(), shots}
    );
  }

auto reserve_response = post(IBM_API_URL, "/runtime/jobs", body.dump(), headers);

}

void IBMAccelerator::cancel() {
  xacc::info("Attempting to cancel IBM job " + currentJobId);
  if (!hub.empty() && jobIsRunning && !currentJobId.empty()) {
    xacc::info("Canceling IBM Job " + currentJobId);
    std::map<std::string, std::string> headers{
        {"Content-Type", "application/x-www-form-urlencoded"},
        {"Connection", "keep-alive"},
        {"Content-Length", "0"}};
    auto path = IBM_CREDENTIALS_PATH + "/Jobs/" + currentJobId +
                "/cancel?access_token=" + currentApiToken;
    auto response = post(IBM_API_URL, path, std::string(), headers);
    xacc::info("Cancel Response: " + response);
    jobIsRunning = false;
    currentJobId = "";
  }
}

std::vector<std::pair<int, int>> IBMAccelerator::getConnectivity() {

  if (!xacc::container::contains(availableBackends, backend)) {
    xacc::error(backend + " is not available.");
  }

  const std::string path("/runtime/backends/" + backend + "/configuration");
  auto backend_json = json::parse(get(IBM_API_URL, path, headers, {{"name", backend}}));
  std::vector<std::pair<int, int>> graph;
  graph = backend_json["coupling_map"].get<std::vector<std::pair<int, int>>>();
  return graph;
}

void IBMAccelerator::searchAPIKey(std::string &key, std::string &hub,
                                  std::string &group, std::string &project) {

  // Search for the API Key in $HOME/.ibm_config,
  // $IBM_CONFIG, or in the command line argument --ibm-api-key
  std::string ibmConfig(std::string(getenv("HOME")) + "/.ibm_config");
  if (xacc::fileExists(ibmConfig)) {
    findApiKeyInFile(key, hub, group, project, ibmConfig);
  } else {
    xacc::error(
        "Cannot find IBM Config file with credentials (~/.ibm_config).");
  }
}

void IBMAccelerator::findApiKeyInFile(std::string &apiKey, std::string &hub,
                                      std::string &group, std::string &project,
                                      const std::string &path) {
  std::ifstream stream(path);
  std::string contents((std::istreambuf_iterator<char>(stream)),
                       std::istreambuf_iterator<char>());

  std::vector<std::string> lines;
  lines = xacc::split(contents, '\n');
  for (auto l : lines) {
    if (l.find("key") != std::string::npos) {
      std::vector<std::string> split = xacc::split(l, ':');
      auto key = split[1];
      xacc::trim(key);
      apiKey = key;
    } else if (l.find("hub") != std::string::npos) {
      std::vector<std::string> split;
      split = xacc::split(l, ':');
      auto _hub = split[1];
      xacc::trim(_hub);
      hub = _hub;
    } else if (l.find("group") != std::string::npos) {
      std::vector<std::string> split;
      split = xacc::split(l, ':');
      auto _group = split[1];
      xacc::trim(_group);
      group = _group;
    } else if (l.find("project") != std::string::npos) {
      std::vector<std::string> split;
      split = xacc::split(l, ':');
      auto _project = split[1];
      xacc::trim(_project);
      project = _project;
    }
  }
}

HeterogeneousMap IBMAccelerator::getProperties() {
  HeterogeneousMap m;

  if (backendProperties.count(backend)) {
    auto props = backendProperties[backend];

    m.insert("total-json", props.dump());
    m.insert("config-json", chosenBackend.dump());
    m.insert("defaults-json", defaults_response);

    auto qubit_props = props["qubits"];

    std::vector<double> p01s, p10s;

    for (auto &qp : qubit_props) {
      for (auto &q : qp) {
        if (q["name"].get<std::string>() == "prob_meas0_prep1") {
          p01s.push_back(q["value"].get<double>());
        } else if (q["name"].get<std::string>() == "prob_meas1_prep0") {
          p10s.push_back(q["value"].get<double>());
        }
      }
    }
    assert(std::all_of(p01s.begin(), p01s.end(),
                       [](double x) { return x >= 0.0 && x < 1.0; }));
    assert(std::all_of(p10s.begin(), p10s.end(),
                       [](double x) { return x >= 0.0 && x < 1.0; }));

    m.insert("p01s", p01s);
    m.insert("p10s", p10s);
    m.insert("multi_meas_enabled", multi_meas_enabled);
  }

  return m;
} // namespace quantum


const std::string RestClient::post(const std::string &remoteUrl,
                                   const std::string &path,
                                   const std::string &postStr,
                                   std::map<std::string, std::string> headers,
                                   const std::string &queryParams) {
#ifndef REMOTE_DISABLED
  if (headers.empty()) {
    headers.insert(std::make_pair("Content-type", "application/json"));
    headers.insert(std::make_pair("Connection", "keep-alive"));
    headers.insert(std::make_pair("Accept", "*/*"));
  }

  cpr::Header cprHeaders;
  for (auto &kv : headers) {
    cprHeaders.insert({kv.first, kv.second});
  }

  if (verbose)
    xacc::info("Posting to " + remoteUrl + path + ", with data " + postStr);

  cpr::Response r;
  if (queryParams.empty()) {
    r = cpr::Post(cpr::Url{remoteUrl + path}, cpr::Body(postStr), cprHeaders,
                     cpr::VerifySsl(false));
  } else {
    json jsonQueryParams = json::parse(queryParams);
    cpr::Parameters cprQueryParams{};
    for (const auto& [key, value] : jsonQueryParams.items()) {
        cprQueryParams.Add({key, value});
    }
        r = cpr::Post(cpr::Url{remoteUrl + path}, cpr::Body(postStr), cprHeaders,
                     cprQueryParams, cpr::VerifySsl(false));

  }

  if (r.status_code != 200)
    throw std::runtime_error("HTTP POST Error - status code " +
                             std::to_string(r.status_code) + ": " +
                             r.error.message + ": " + r.text);

  return r.text;
#else
  return "";
#endif
}

void RestClient::put(const std::string &remoteUrl, const std::string &putStr,
                     std::map<std::string, std::string> headers) {
#ifndef REMOTE_DISABLED
  if (headers.empty()) {
    headers.insert(std::make_pair("Content-type", "application/json"));
    headers.insert(std::make_pair("Connection", "keep-alive"));
    headers.insert(std::make_pair("Accept", "*/*"));
  }

  cpr::Header cprHeaders;
  for (auto &kv : headers) {
    cprHeaders.insert({kv.first, kv.second});
  }

  if (verbose)
    xacc::info("PUT to " + remoteUrl + " with data " + putStr);
  auto r = cpr::Put(cpr::Url{remoteUrl}, cpr::Body(putStr), cprHeaders,
                    cpr::VerifySsl(false));

  if (r.status_code != 200)
    throw std::runtime_error("HTTP POST Error - status code " +
                             std::to_string(r.status_code) + ": " +
                             r.error.message + ": " + r.text);
#endif
  return;
}
const std::string
RestClient::get(const std::string &remoteUrl, const std::string &path,
                std::map<std::string, std::string> headers,
                std::map<std::string, std::string> extraParams) {
#ifndef REMOTE_DISABLED
  if (headers.empty()) {
    headers.insert(std::make_pair("Content-type", "application/json"));
    headers.insert(std::make_pair("Connection", "keep-alive"));
    headers.insert(std::make_pair("Accept", "*/*"));
  }

  cpr::Header cprHeaders;
  for (auto &kv : headers) {
    cprHeaders.insert({kv.first, kv.second});
  }

  cpr::Parameters cprParams;
  for (auto &kv : extraParams) {
    cprParams.Add({kv.first, kv.second});
  }

  if (verbose)
    xacc::info("GET at " + remoteUrl + path);
  auto r = cpr::Get(cpr::Url{remoteUrl + path}, cprHeaders, cprParams,
                    cpr::VerifySsl(false));

  if (r.status_code != 200)
    throw std::runtime_error("HTTP GET Error - status code " +
                             std::to_string(r.status_code) + ": " +
                             r.error.message + ": " + r.text);

  return r.text;
#else
  return "";
#endif
}

std::string IBMAccelerator::post(const std::string &_url,
                                 const std::string &path,
                                 const std::string &postStr,
                                 std::map<std::string, std::string> headers,
                                 const std::string& queryParams) {
  std::string postResponse;
  int retries = 10;
  std::exception ex;
  bool succeeded = false;

  // Execute HTTP Post
  do {
    try {
      postResponse = restClient->post(_url, path, postStr, headers, queryParams);
      succeeded = true;
      break;
    } catch (std::exception &e) {
      ex = e;
      xacc::info("Remote Accelerator " + name() +
                 " caught exception while calling restClient->post() "
                 "- " +
                 std::string(e.what()));
      retries--;
      if (retries > 0) {
        xacc::info("Retrying HTTP Post.");
      }
    }
  } while (retries > 0);

  if (!succeeded) {
    cancel();
    xacc::error("Remote Accelerator " + name() +
                " failed HTTP Post for Job Response - " +
                std::string(ex.what()));
  }

  return postResponse;
}

void IBMAccelerator::put(const std::string &_url, const std::string &postStr,
                         std::map<std::string, std::string> headers) {
  int retries = 10;
  std::exception ex;
  bool succeeded = false;

  // Execute HTTP Post
  do {
    try {
      restClient->put(_url, postStr, headers);
      succeeded = true;
      break;
    } catch (std::exception &e) {
      ex = e;
      xacc::info("Remote Accelerator " + name() +
                 " caught exception while calling restClient->put() "
                 "- " +
                 std::string(e.what()));
      retries--;
      if (retries > 0) {
        xacc::info("Retrying HTTP Put.");
      }
    }
  } while (retries > 0);

  if (!succeeded) {
    cancel();
    xacc::error("Remote Accelerator " + name() +
                " failed HTTP Put for Job Response - " +
                std::string(ex.what()));
  }

  return;
}

// TODO
std::string
IBMAccelerator::getNativeCode(std::shared_ptr<CompositeInstruction> program,
                              const HeterogeneousMap &config) {
  return "";
}

std::string
IBMAccelerator::get(const std::string &_url, const std::string &path,
                    std::map<std::string, std::string> headers,
                    std::map<std::string, std::string> extraParams) {
  std::string getResponse;
  int retries = 10;
  std::exception ex;
  bool succeeded = false;
  // Execute HTTP Get
  do {
    try {
      getResponse = restClient->get(_url, path, headers, extraParams);
      succeeded = true;
      break;
    } catch (std::exception &e) {
      ex = e;
      xacc::info("Remote Accelerator " + name() +
                 " caught exception while calling restClient->get() "
                 "- " +
                 std::string(e.what()));
      // s1.find(s2) != std::string::npos) {
      if (std::string(e.what()).find("Caught CTRL-C") != std::string::npos) {
        cancel();
        xacc::error(std::string(e.what()));
      }
      retries--;
      if (retries > 0) {
        xacc::info("Retrying HTTP Get.");
      }
    }
  } while (retries > 0);

  if (!succeeded) {
    cancel();
    xacc::error("Remote Accelerator " + name() +
                " failed HTTP Get for Job Response - " +
                std::string(ex.what()));
  }

  return getResponse;
}

} // namespace quantum
} // namespace xacc