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

// This provides the implementation for the Qubit-Wise Commutativity grouping

#ifndef XACC_OBSERVABLE_GROUPING_QWC_HPP_
#define XACC_OBSERVABLE_GROUPING_QWC_HPP_

#include "ObservableGrouping.hpp"
#include "xacc.hpp"
#include "PauliOperator.hpp"

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {

class QubitWiseCommutativity : public ObservableGrouping {

protected:
  // we may want to move this as a pure virtual method 
  // in ObservableGrouping as each implementation has 
  // a commutativity condition, but will leave this here for now
  bool checkIfCommute(const std::map<int, std::string> op1,
                      const std::map<int, std::string> op2) const;

public:
  std::vector<std::shared_ptr<Observable>>
  group(std::shared_ptr<Observable> observable) override;

  std::vector<std::shared_ptr<CompositeInstruction>>
  observeGroups(std::vector<std::shared_ptr<Observable>> groups) override;

  const std::string name() const override { return "qwc"; }
  const std::string description() const override { return ""; }
};

} // namespace xacc
#endif