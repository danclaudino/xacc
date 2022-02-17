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

#ifndef XACC_ACCELERATOR_BUFFER_DECORATOR_HPP_
#define XACC_ACCELERATOR_BUFFER_DECORATOR_HPP_

#include "AcceleratorBuffer.hpp"
#include "Identifiable.hpp"

namespace xacc {

// we make this inherit from Identifiable in case there may be more such codes
class AcceleratorBufferDecorator : public AcceleratorBuffer, public Identifiable {
protected:
  std::shared_ptr<AcceleratorBuffer> decoratedAcceleratorBuffer;


public:
  AcceleratorBufferDecorator() {}
  AcceleratorBufferDecorator(std::shared_ptr<AcceleratorBuffer> a)
      : decoratedAcceleratorBuffer(a) {}

  void setDecorated(std::shared_ptr<AcceleratorBuffer> a) {
    decoratedAcceleratorBuffer = a;
  }

  void addExtraInfo(const std::string infoName, ExtraInfo i) {
    decoratedAcceleratorBuffer->addExtraInfo(infoName, i);
    return;
  }

  bool hasExtraInfoKey(const std::string infoName) {
  return decoratedAcceleratorBuffer->getInformation().count(infoName);
  }

  void appendChild(const std::string name,
                   std::shared_ptr<AcceleratorBuffer> buffer) {
  decoratedAcceleratorBuffer->appendChild(name, buffer);
  }

  const int size() { return decoratedAcceleratorBuffer->size(); }

  std::shared_ptr<AcceleratorBuffer> getDecorated() { return decoratedAcceleratorBuffer; }
};

} // namespace xacc

#endif