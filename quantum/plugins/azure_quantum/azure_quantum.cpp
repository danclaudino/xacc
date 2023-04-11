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
#include "azure_quantum.hpp"

#include <cpr/cpr.h>
#include <cstdlib>
#include <spdlog/fmt/fmt.h>

#include <cctype>
#include <ctime>
#include <fstream>
#include <thread>

#include "Utils.hpp"
#include "json.hpp"
#include "xacc_service.hpp"
#include "xacc_plugin.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace xacc {
namespace quantum {

void AzureQuantumAccelerator::initialize(const HeterogeneousMap &params) {
  if (!initialized) {
    const std::string azureConfigFilename(std::string(getenv("HOME")) +
                                          "/.azure_config");

    if (!xacc::fileExists(azureConfigFilename)) {
      xacc::error("Could not find Azure Quantum configuration file " +
                  azureConfigFilename);
    }

    if (params.keyExists<int>("shots")) {
      shots = params.get<int>("shots");
    }

    if (params.stringExists("job-name")) {
      job_name = params.getString("job-name");
    }

    if (params.stringExists("retrieve-job-id")) {
      retrieve_job_id = params.getString("job-id");
    }
    auto configInfo = azureUtils::getConfigInfo(azureConfigFilename);
    baseUrl = configInfo.first;
    accessToken = configInfo.second;

    initialized = true;
  }
}

void AzureQuantumAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> circuit) {
      /*
  if (backend.empty()) {
    xacc::error(
        "Please specify a backend in your getAccelerator() call.");
  }
  */

    // At a high level we need to:
    // 1. Authenticate using AAD
    // This is required to get an access token (in the .azure_config file)
    // "qcor -set-credentials azure"
        // std::cout << "baseUrl: " << baseUrl << std::endl;
        // std::cout << "token: " << accessToken << std::endl;
        cpr::Header cprHeaders;
        accessToken = "eyJhbGciOiJSU0EtT0FFUCIsImVuYyI6IkExMjhDQkMtSFMyNTYiLCJ4NXQiOiJZUEVTNE1zc090NjJYLWxnX3dMcHJyc2c4RG8iLCJ6aXAiOiJERUYifQ.MdbuDQM3iDwDvDF61w_jlKm-nk00zlL1l2ap0xUCM_k1Sc9K4eF1rpgQHLRdGpaPvfLq0si97VV2fqDsTtfhvfjVtC-iyL5qYFsSxmffFvLsBy5xHyOrZ39G6d-c4Sa3vMLUGw8oKmCdUWu-7jyHdjUtA6d4nflIHTUB2zR2AzQuRcQe16eQvaLK0S4O7kEdtNrErYQhVM1wA1tfXUwEAaqTGHGZiBn3PSZcJKnrgy9fudUjHqHRKsXoB9LtPHl13elMeem3M-poeSYa1idTK8Indl82FT71yAjMIkEN1N1Rr8C43lu_Z3G57O6KxhLc17ekUXw0ZDshH9VVtDbPaQ.xy4QchCACSPw4tuamP8spQ.UXCgBZP0iWPCZMjUCf1qQx2Q8h7xcIKryEjZNR8aBXQTpl62gUyUQ3n_fnFTylbocUzMxTFVUVYtWyJoZVOuv3P3jGC3NeK4LBARLct_Xpan9CZrT-zD5ElHq5uUqf6PvGqQiaACRRZ1xn4KjnG5E4MfQ9Vd0dvOCerqm6rW_zkYTHZguAhCF-isZSAEec0GXLlsNdgfQ1KOmV5IKt-JyYrX6euNi2_yb2ddDlLBGHlN0KyImPqWyIJ_5IMVNGN2ndK7D9WU0x-06h0YXHsXU5TwtryohC_2ye3W66JfNbXD5XL6IwCG9pKB83-Qq_e_A4V-NqkbheV1mVIQ_wEx44v7kidHXukUUsJEBGWovTqgro6PeR8PIR_RxZ7eVzeoQXZl-dx9kxJ5aang2ExGSJT1cu5xiDEVAEqCq9aVEpCTTNIK7fuYd2du-1HyPOGVjiMTCFr3-YeypsBxqUanKDMMpTeM7Bw2pa1nofNwsNgLN6QCML630sA7eHzBJHMw3X5uqxqHocsY-OXU9cVpTx1Nn_0S1feawHPQoJjy7FPKGyv5aQnOu2uhbKk5HLv95Zj0s5rNXhK_vpHS5veF4S1MxK4FgCOtC-SwwuWrOaCCLKlj-cvGeoVzN4FB6mTAT8Ol2V-8pTdhpjZagB_jJIwVyMEXddIBSj4Tg1KON_GB09cA0P0ZRysbrK3W_WDMr_hKSSpC5E4erOa0oWPykXiIuJFKFbGaUkDFZinZr3I1jr5A0UkR3xndUfX0_Q34jHuDiW2h0UM9QXfm13gk66lLCDFc3Or_TntF5M56PX1-YUWX0EB4WyG-e3ZJFB2nJF_ZDLBVZQ7a653tqMuy1KzyIauF9rUUjC0P7RMGWRiXMly4E0ft8g15ZVTEwP8qQYGQ-ubw1UfG2RTIqQAwvdEMTo2FXcCHqNkoCBAGrxccl2Mm12oJL76RTt39CV8k2hmcl18s2NblqxI1yEtwMdQT_yPLjEqHGCAKkBWtPOScbHJ1_7ehTAjGN1hxuff1hQGDMseJsqhKAxqDuyC31fng9t2ESJxwN9PFLONinbKnSwtSDXOfo8N5KBxfnuwvviukGqrxUOsGy1VKSxhsIm0lQYw6lhwJxesOHSSGxOc.6yZyobBJ1efBE9W-yXh4bg";
        cprHeaders.insert({"Content-type", "application/json"});
        cprHeaders.insert({"Connection", "keep-alive"});
        cprHeaders.insert({"Accept", "application/json"});
        // Please replace the User-Agent value with a descriptive name for your client, e.g. "qcor-cpp".
        cprHeaders.insert({"User-Agent", "azure-quantum-cpp-demo-v0.1"});
        cprHeaders.insert({"Authorization", "Bearer " + accessToken});

        // Testing: query provider status
        const std::string path = "/providerStatus";
        cpr::Parameters cprParams;
        auto sasResponse = cpr::Get(cpr::Url{baseUrl + path}, cprHeaders, cprParams,
                                    cpr::VerifySsl(false));
        std::cout << "\nResponse:" << sasResponse.text << "\n";
  //std::cout << accessToken <<  " HERE\n";
  exit(0);
        // 2. Upload the QIR definition to an Azure Storage
        // see
        // https://docs.microsoft.com/en-us/rest/api/storageservices/put-blob-from-url
        // Get the SAS token for the Azure Storage storage account associated with
        // the quantum workspace. The SAS URL can be used to upload job input and/or
        // download job output.
        const std::string storagePath = baseUrl + "/storage/sasUri";
        const auto jobId = azureUtils::randomIdStr();
        const std::string jobName = "qcor-" + jobId;

        // POST
        nlohmann::json j;
        j["containerName"] = jobName;
        auto r = cpr::Post(cpr::Url{storagePath}, cprHeaders, cpr::Body{j.dump()},
                           cpr::VerifySsl(false));
        std::cout << "SAS Response:" << r.text << "\n";
        nlohmann::json sasJson = nlohmann::json::parse(r.text);
        const std::string sasUri = sasJson["sasUri"].get<std::string>();

        // PUT the JOB data to the Azure Storage
        std::ifstream ifs("Sample.bc");
        std::string content;
        content.assign( (std::istreambuf_iterator<char>(ifs) ),
                        (std::istreambuf_iterator<char>()    ) );

        const auto [eTag, requestId] =
            createContainerBlob(azureUtils::Url{sasUri});
        std::cout << "Etag:" << eTag << "\n";
        std::stringstream stream, gzipped;
        stream << content;
        static constexpr const char *blobName = "inputData";

        // GZIP before uploading
        boost::iostreams::filtering_streambuf< boost::iostreams::input> in;
        in.push(boost::iostreams::gzip_compressor());
        in.push(stream);
        boost::iostreams::copy(in, gzipped);

        // Upload data and get blob URI
        const auto blobUri =
            uploadContainerBlob(azureUtils::Url{sasUri}, blobName, gzipped);

        // For testing: use an existsing blob URI
        // const auto blobUri = "https://e2etests.blob.core.windows.net/job-94175976-8a3e-11ec-b6ac-00155d6d8f85/inputData";

        std::cout << "Blob access url:" << blobUri << std::endl;

        // 3. Create the job metadata and submit a job request to Azure.
        // https://docs.microsoft.com/en-us/rest/api/azurequantum/dataplane/jobs
        
        const std::string path1 = "/jobs/" + jobId;
        nlohmann::json job;
        nlohmann::json inputParams;
        inputParams["entryPoint"] = "Sample__HelloQ";

        job["id"] = jobId;
        job["name"] = "QIR test";
        job["containerUri"] = sasUri;
        job["inputDataUri"] = blobUri;
        job["inputDataFormat"] = "qir.v1/full-profile";
        job["outputDataFormat"] = "microsoft.qio-results.v2";
        job["inputParams"] = inputParams;
        job["providerId"] = "Microsoft.Simulator";
        job["target"] = "microsoft.simulator.fullstate";
        job["outpuDataFormat"] = "microsoft.qio-results.v2";
        auto sasResponse1 = cpr::Put(
            cpr::Url{baseUrl + path1},
            cprHeaders,
            cpr::Body{job.dump()},
            cpr::VerifySsl(false)
        );
        std::cout << "Job details:\n" << sasResponse1.text << "\n";

        // 4.a Get job status
        auto sasResponse2 = cpr::Get(
            cpr::Url{baseUrl + path1},
            cprHeaders,
            cprParams,
            cpr::VerifySsl(false)
        );
        nlohmann::json jobDetails = nlohmann::json::parse(sasResponse2.text);
        std::cout << "\nJob details:\n" << sasResponse2.text << "\n";

        // 4.b Wait for job to complete
        if (jobDetails["status"] == "Waiting") {
            std::thread t{azureUtils::wait_until_done};
            t.join();
        }

        // 5. Fetch results from container
        static constexpr const char *outputBlobName = "rawOutputData";
        auto output = downloadContainerBlob(azureUtils::Url{sasUri}, outputBlobName);
        std::cout << "\nResults:\n" << output;

  return;
}

void AzureQuantumAccelerator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>>
        compositeInstructions) {
  for (auto &f : compositeInstructions) {
    auto tmpBuffer =
        std::make_shared<xacc::AcceleratorBuffer>(f->name(), buffer->size());
    execute(tmpBuffer, f);
    buffer->appendChild(f->name(), tmpBuffer);
  }
}

/*
std::string
AzureQuantumAccelerator::post(const std::string &_url, const std::string &path,
                              const std::string &postStr,
                              std::map<std::string, std::string> headers) {
  std::string postResponse;
  int retries = 10;
  std::exception ex;
  bool succeeded = false;

  // Execute HTTP Post
  do {
    try {
      postResponse = restClient->post(_url, path, postStr, headers);
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

std::string
AzureQuantumAccelerator::get(const std::string &_url, const std::string &path,
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


std::map<std::string, std::string>
AzureQuantumAccelerator::generateRequestHeader(const std::string user_agent) const {
  std::map<std::string, std::string> headers{
      {"Authorization", accessToken},
      {"User-Agent", user_agent},
      {"Content-Type", "application/json"},
      {"Connection", "keep-alive"},
      {"Accept", "application/json"}};
  return headers;
}
*/

} // namespace quantum
} // namespace xacc

REGISTER_ACCELERATOR(xacc::quantum::AzureQuantumAccelerator)