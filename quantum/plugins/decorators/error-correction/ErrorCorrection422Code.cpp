#include "ErrorCorrection422Code.hpp"
#include <numeric>
#include "xacc.hpp"

namespace xacc {
namespace quantum {

void ErrorCorrection422Code::appendMeasurement(const std::string &measurement) {
  
  // this assumes the first qubit is the post selecting ancilla
  // and the last determines theta

  int ancilla0;
  /*
  if(this->hasExtraInfoKey("bit-order")) {
  if (this->getInformation("bit-order") == 1) {
    // MSB
    ancilla0 = std::stoi(&measurement.front());
  } else {
    ancilla0 = std::stoi(&measurement.back());
  }
  }
*/

  ancilla0 = checkAncillaMeasurement(measurement, 0);

  // value 1 in the first ancilla are discarded
  if (!ancilla0) {
    postSelectedBitStringToCounts[measurement]++;
  }

  return;

}

void ErrorCorrection422Code::appendMeasurement(const std::string measurement, const int count) {
  
  // this assumes the first qubit is the post selecting ancilla
  // and the last determines theta

  int ancilla0;
  /*
  if(this->hasExtraInfoKey("bit-order")) {
  if (this->getInformation("bit-order") == 1) {
    // MSB
    ancilla0 = std::stoi(&measurement.front());
  } else {
    ancilla0 = std::stoi(&measurement.back());
  }
  }
  */
  ancilla0 = checkAncillaMeasurement(measurement, 0);

  // value 1 in the first ancilla are discarded
  if (!ancilla0) {
    postSelectedBitStringToCounts[measurement] = count;
  }

  return;

}

double
ErrorCorrection422Code::computeMeasurementProbability(const std::string &bitStr) {
  return (double)postSelectedBitStringToCounts[bitStr] /
         std::accumulate(
             postSelectedBitStringToCounts.begin(), postSelectedBitStringToCounts.end(), 0,
             [](int value, const std::map<std::string, int>::value_type &p) {
               return value + p.second;
             });
}

const double ErrorCorrection422Code::getExpectationValueZ() {
  double aver = 0.0;
  std::pair<double, double> avers = {0.0, 0.0};
  auto has_even_parity = [](unsigned int x) -> int {
    unsigned int count = 0, i, b = 1;
    for (i = 0; i < 32; i++) {
      if (x & (b << i)) {
        count++;
      }
    }
    if ((count % 2)) {
      return 0;
    }
    return 1;
  };

  // check the second ancilla

  if (this->hasExtraInfoKey("exp-val-z")) {
    aver = mpark::get<double>(getInformation("exp-val-z"));
  } else {
    if (postSelectedBitStringToCounts.empty()) {
      xacc::error("called getExpectationValueZ() on an AcceleratorBuffer with "
                  "no measurements!");
      return 0;
    }

    for (auto &kv : postSelectedBitStringToCounts) {
      int i = std::stoi(kv.first, nullptr, 2);
      auto ancilla2 = checkAncillaMeasurement(kv.first, 1);
      auto par = has_even_parity(i);
      auto p = computeMeasurementProbability(kv.first);
      if (!par) {
        p = -p;
      }
      if (ancilla2) {
        avers.second += p;
      } else {
        avers.first += p;
      }
      aver += p;
    }
  }
  return aver;
}

} // namespace quantum
} // namespace quantum