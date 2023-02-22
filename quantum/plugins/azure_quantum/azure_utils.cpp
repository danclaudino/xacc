#include "azure_utils.hpp"
#include <algorithm>
#include <cctype>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <cassert>
#include <cpr/cpr.h>
#include <random>
#include <thread>

namespace {
std::string FormatEncodedUrlQueryParameters(
    std::map<std::string, std::string> const &encodedQueryParameters) {
  {
    std::string queryStr;
    if (!encodedQueryParameters.empty()) {
      auto separator = '?';
      for (const auto &q : encodedQueryParameters) {
        queryStr += separator + q.first + '=' + q.second;
        separator = '&';
      }
    }

    return queryStr;
  }
}
} // namespace

// namespace qcor {
namespace azureUtils {

std::string Url::GetUrlWithoutQuery(bool relative) const {
  std::string url;

  if (!relative) {
    if (!m_scheme.empty()) {
      url += m_scheme + "://";
    }
    url += m_host;
    if (m_port != 0) {
      url += ":" + std::to_string(m_port);
    }
  }

  if (!m_encodedPath.empty()) {
    if (!relative) {
      url += "/";
    }

    url += m_encodedPath;
  }

  return url;
}

void Url::AppendQueryParameters(const std::string &query) {
  std::string::const_iterator cur = query.begin();
  if (cur != query.end() && *cur == '?') {
    ++cur;
  }

  while (cur != query.end()) {
    auto key_end = std::find(cur, query.end(), '=');
    std::string query_key = std::string(cur, key_end);

    cur = key_end;
    if (cur != query.end()) {
      ++cur;
    }

    auto value_end = std::find(cur, query.end(), '&');
    std::string query_value = std::string(cur, value_end);

    cur = value_end;
    if (cur != query.end()) {
      ++cur;
    }
    m_encodedQueryParameters[std::move(query_key)] = std::move(query_value);
  }
}

Url::Url(const std::string &url) {
  std::string::const_iterator pos = url.begin();
  const std::string schemeEnd = "://";
  auto schemeIter = url.find(schemeEnd);
  if (schemeIter != std::string::npos) {
    std::transform(url.begin(), url.begin() + schemeIter,
                   std::back_inserter(m_scheme), ::tolower);

    pos = url.begin() + schemeIter + schemeEnd.length();
  }

  auto hostIter = std::find_if(
      pos, url.end(), [](char c) { return c == '/' || c == '?' || c == ':'; });
  m_host = std::string(pos, hostIter);
  pos = hostIter;

  if (pos != url.end() && *pos == ':') {
    auto port_ite = std::find_if_not(pos + 1, url.end(), [](char c) {
      return std::isdigit(static_cast<unsigned char>(c));
    });
    auto portNumber = std::stoi(std::string(pos + 1, port_ite));

    // stoi will throw out_of_range when `int` is overflow, but we need to throw
    // if uint16 is overflow
    auto maxPortNumberSupported = std::numeric_limits<uint16_t>::max();
    if (portNumber > maxPortNumberSupported) {
      throw std::out_of_range(
          "The port number is out of range. The max supported number is " +
          std::to_string(maxPortNumberSupported) + ".");
    }
    // cast is safe because the overflow was detected before
    m_port = static_cast<uint16_t>(portNumber);
    pos = port_ite;
  }

  if (pos != url.end() && (*pos != '/') && (*pos != '?')) {
    // only char `\` or `?` is valid after the port (or the end of the URL). Any
    // other char is an invalid input
    throw std::invalid_argument("The port number contains invalid characters.");
  }

  if (pos != url.end() && (*pos == '/')) {
    auto pathIter = std::find(pos + 1, url.end(), '?');
    m_encodedPath = std::string(pos + 1, pathIter);
    pos = pathIter;
  }

  if (pos != url.end() && *pos == '?') {
    auto queryIter = std::find(pos + 1, url.end(), '#');
    AppendQueryParameters(std::string(pos + 1, queryIter));
    pos = queryIter;
  }
}

std::string Url::GetRelativeUrl() const {
  return GetUrlWithoutQuery(true) +
         FormatEncodedUrlQueryParameters(m_encodedQueryParameters);
}
std::string Url::GetAbsoluteUrl() const {
  return GetUrlWithoutQuery(false) +
         FormatEncodedUrlQueryParameters(m_encodedQueryParameters);
}

std::pair<std::string, std::string> createContainerBlob(const Url &url) {
  auto request = url;
  request.AppendQueryParameter("restype", "container");
  cpr::Header cprHeaders;
  cprHeaders.insert({"Content-Length", "0"});
  cprHeaders.insert({"x-ms-version", "2020-02-10"});
  auto response = cpr::Put(cpr::Url{request.GetAbsoluteUrl()}, cprHeaders,
                           cpr::VerifySsl(false));
  if (response.status_code != 201) {
    throw std::runtime_error("Failed to create container blob");
  }
  return std::make_pair(response.header["ETag"],
                        response.header["x-ms-request-id"]);
}

std::string uploadContainerBlob(const Url &url, const std::string &blobName,
                                std::stringstream &stream) {
  auto request = url;
  request.AppendPath(blobName);
  cpr::Header cprHeaders;
  const auto requestBody = stream.str();
  cprHeaders.insert({"Content-Length", std::to_string(requestBody.size())});
  cprHeaders.insert({"x-ms-version", "2020-02-10"});
  cprHeaders.insert({"Accept-Encoding", "gzip, deflate"});
  cprHeaders.insert({"Connection", "keep-alive"});
  cprHeaders.insert({"x-ms-blob-type", "BlockBlob"});
  cprHeaders.insert({"Accept", "application/xml"});
  cprHeaders.insert({"x-ms-blob-content-encoding", "gzip"});
  cprHeaders.insert({"x-ms-blob-content-type", "application/json"});
  cprHeaders.insert({"Content-Type", "application/octet-stream"});
  auto response = cpr::Put(cpr::Url{request.GetAbsoluteUrl()}, cprHeaders,
                           cpr::Body{requestBody}, cpr::VerifySsl(false));
  if (response.status_code != 201) {
    throw std::runtime_error("Failed to upload container blob");
  }

  return request.GetUrlWithoutQuery(false);
}

std::string downloadContainerBlob(const Url &url, const std::string &blobName) {
  auto request = url;
  request.AppendPath(blobName);
  cpr::Header cprHeaders;
  cpr::Parameters cprParams;
  cprHeaders.insert({"x-ms-version", "2020-02-10"});
  cprHeaders.insert({"Accept-Encoding", "gzip, deflate, br"});
  cprHeaders.insert({"Connection", "keep-alive"});
  cprHeaders.insert({"x-ms-blob-type", "BlockBlob"});
  cprHeaders.insert({"Accept", "application/xml"});
  cprHeaders.insert({"x-ms-blob-content-encoding", "gzip"});
  cprHeaders.insert({"x-ms-blob-content-type", "application/json"});
  cprHeaders.insert({"Content-Type", "application/octet-stream"});
  auto response = cpr::Get(cpr::Url{request.GetAbsoluteUrl()}, cprHeaders,
                           cprParams, cpr::VerifySsl(false));
  return response.text;
}

std::string getBaseUrl(const std::string &region,
                       const std::string &subscriptionId,
                       const std::string &resourceGroup,
                       const std::string &workspaceName) {
  const auto formatRegion = [](const std::string &rawRegionString) {
    std::string regionString = rawRegionString;
    // e.g., East US -> eastus
    regionString.erase(
        std::remove_if(regionString.begin(), regionString.end(), isspace),
        regionString.end());
    std::transform(regionString.begin(), regionString.end(),
                   regionString.begin(), ::tolower);
    return regionString;
  };
  return "https://" + formatRegion(region) +
         ".quantum.azure.com/v1.0/subscriptions/" + subscriptionId +
         "/resourceGroups/" + resourceGroup +
         "/providers/Microsoft.Quantum/workspaces/" + workspaceName;
}

std::pair<std::string, std::string>
getConfigInfo(const std::string &configFilePath) {

  std::ifstream configFile(configFilePath);
  const std::string contents((std::istreambuf_iterator<char>(configFile)),
                             std::istreambuf_iterator<char>());
  const auto lines = xacc::split(contents, '\n');
  std::string token, location, subscriptionId, resourceGroup, workspaceName;

  for (const auto &line : lines) {
    if (line.find("key") != std::string::npos) {
      xacc::trim(xacc::split(line, ':')[1]);
      token = line;
    }
    if (line.find("Location") != std::string::npos) {
      xacc::trim(xacc::split(line, ':')[1]);
      location = line;
    }
    if (line.find("expires") != std::string::npos) {
      xacc::trim(xacc::split(line, ':')[1]);
      const auto expiration = std::stoll(line);
      if (std::time(nullptr) > expiration) {
        throw std::runtime_error("Token expired");
      }
    }
    if (line.find("Resource ID") != std::string::npos) {
      xacc::trim(xacc::split(line, ':')[1]);
      const auto resourcePath = line;
      const auto components = xacc::split(resourcePath, '/');
      for (size_t i = 0; i < components.size(); ++i) {
        if (components[i] == "subscriptions") {
          subscriptionId = components[i + 1];
        }
        if (components[i] == "resourceGroups") {
          resourceGroup = components[i + 1];
        }
        if (components[i] == "Workspaces") {
          workspaceName = components[i + 1];
        }
      }
    }
  }
  
  assert(!token.empty());
  assert(!subscriptionId.empty());
  assert(!resourceGroup.empty());
  assert(!workspaceName.empty());
  assert(!location.empty());
  return std::make_pair(
      getBaseUrl(location, subscriptionId, resourceGroup, workspaceName),
      token);
}

std::string randomIdStr() {
  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::uniform_int_distribution<int> dist(0, 15);

  const char *v = "0123456789abcdef";
  const bool dash[] = {0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0};

  std::string res;
  for (int i = 0; i < 16; ++i) {
    if (dash[i]) {
      res += "-";
    }
    res += v[dist(rng)];
    res += v[dist(rng)];
  }
  return res;
}

void wait(int seconds)
{
  std::this_thread::sleep_for(std::chrono::seconds{seconds});
}

void wait_until_done()
{
  for (int i = 0; i < 20; ++i)
  {
    wait(1);
    std::cout << '.';
  }
}

} // namespace utils
  //} // namespace qcor