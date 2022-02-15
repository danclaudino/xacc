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

namespace xacc {

class AcceleratorBufferDecorator : public AcceleratorBuffer {
protected:
  std::shared_ptr<AcceleratorBuffer> decoratedAcceleratorBuffer;

public:
  AcceleratorBufferDecorator() {}
  AcceleratorBufferDecorator(std::shared_ptr<AcceleratorBuffer> a)
      : decoratedAcceleratorBuffer(a) {}
  void setDecorated(std::shared_ptr<AcceleratorBuffer> a) {
    decoratedAcceleratorBuffer = a;
  }
  std::shared_ptr<AcceleratorBuffer> getDecorated() { return decoratedAcceleratorBuffer; }
};

} // namespace xacc

#endif