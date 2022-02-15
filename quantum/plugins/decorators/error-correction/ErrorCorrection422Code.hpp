/*******************************************************************************
 * Copyright (c) 2022 UT-Battelle, LLC.
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
#ifndef XACC_ERROR_CORRECTION_422_DECORATOR_HPP_
#define XACC_ERROR_CORRECTION_422_DECORATOR_HPP_

#include "AcceleratorBufferDecorator.hpp"
#include "AcceleratorDecorator.hpp"
#include <cstddef>

namespace xacc {
namespace quantum {

class ErrorCorrection422Code : public AcceleratorBufferDecorator {

private:

// if ancilla = 0, then it'spost-selection
// if ancilla = 1, then it's the theta ancilla
bool checkAncillaMeasurement(const std::string bitString, const int ancilla) {

  int bitOrder = 0;
  if(this->hasExtraInfoKey("bit-order")) {
    bitOrder = this->getInformation("bit-order").as<int>();
  }

  if(bitOrder) {
    if (ancilla) return std::stoi(&bitString.back());
    return std::stoi(&bitString.front());
  } else {
    if (ancilla) return std::stoi(&bitString.front());
    return std::stoi(&bitString.back());
  }

};

protected:
  std::map<std::string, int> postSelectedBitStringToCounts;
  // this is the assumed default LSB, always a good idea to provide it
  int bitOrder = 0;

public:

  void appendMeasurement(const std::string &measurement) override;
  void appendMeasurement(const std::string measurement,
                                 const int count) override;
  double computeMeasurementProbability(const std::string &bitStr) override;

  const double getExpectationValueZ() override;
};

} // namespace quantum
} // namespace xacc


#endif