/*******************************************************************************
 * Copyright (c) 2020 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Daniel Claudino - initial API and implementation
 ******************************************************************************/

#ifndef XACC_IR_OBSERVABLE_GROUPING_HPP_
#define XACC_IR_OBSERVABLE_GROUPING_HPP_

#include "Observable.hpp"
#include <vector>

namespace xacc {

class ObservableGrouping : public Identifiable {
public:
  virtual std::vector<std::shared_ptr<Observable>>
  group(std::shared_ptr<Observable> observable) = 0;

  virtual std::vector<std::shared_ptr<CompositeInstruction>>
  observeGroups(std::vector<std::shared_ptr<Observable>>) = 0;
};

} // namespace xacc
#endif