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

#include "qwc.hpp"

#include "PauliOperator.hpp"
#include "xacc.hpp"
#include "xacc_plugin.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include <map>
#include <memory>
#include <vector>
#include "ObservableTransform.hpp"
#include <Eigen/Dense>

using namespace xacc;
using namespace xacc::quantum;

namespace xacc {

std::vector<std::shared_ptr<Observable>>
QubitWiseCommutativity::group(std::shared_ptr<Observable> observable) {

  // first make sure observable is of PauliOperator type
  if (!std::dynamic_pointer_cast<PauliOperator>(observable)) {
    observable =
        xacc::getService<ObservableTransform>("jw")->transform(observable);
  }

  std::vector<std::vector<int>> vertexIdx;
  auto terms = observable->getNonIdentitySubTerms();
  for (int i = 0; i < terms.size(); i++) {

    std::vector<int> v;

    // all of this to get the map<int, string>
    auto term1 = std::dynamic_pointer_cast<PauliOperator>(terms[i])
                     ->getTerms()
                     .begin()
                     ->second.ops();

    for (int j = 0; j < terms.size(); j++) {

      auto term2 = std::dynamic_pointer_cast<PauliOperator>(terms[j])
                       ->getTerms()
                       .begin()
                       ->second.ops();

      v.push_back(checkIfCommute(term1, term2));
    }
    vertexIdx.push_back(v);
  }

  // these are the qwc groups
  std::vector<std::shared_ptr<Observable>> groups;
  // indices of all operators that qw commute
  std::set<std::vector<int>> s(vertexIdx.begin(), vertexIdx.end());

  for (auto i : s) {
    PauliOperator commutingGroup;
    for (int k = 0; k < vertexIdx.size(); k++) {
      if (i == vertexIdx[k]) {
        commutingGroup += *(std::dynamic_pointer_cast<PauliOperator>(terms[k]));
      }
    }
    groups.push_back(std::make_shared<PauliOperator>(commutingGroup));
  }

  return groups;
}

std::vector<std::shared_ptr<CompositeInstruction>>
QubitWiseCommutativity::observeGroups(
    std::vector<std::shared_ptr<Observable>> groups) {

  std::vector<std::shared_ptr<CompositeInstruction>> groupMeasurements;

  for (int group = 0; group < groups.size(); group++) {

    auto provider = xacc::getIRProvider("quantum");
    auto measurements = provider->createComposite("measurements");
    // auto pauliOp = std::dynamic_pointer_cast<PauliOperator>(groups[group]);
    auto obsTerms = groups[group]->getSubTerms();
    auto nOps = obsTerms.size();

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(nOps, nOps);
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(nOps, nOps);

    // Construct Z and X matrices
    // If a term is all Z's, then the X matrix will not be full rank
    // allZ keeps track of where this happens
    std::vector<int> allZ;
    for (int i = 0; i < nOps; i++) {

      int nZeros = 0;
      auto pauliOp = std::dynamic_pointer_cast<PauliOperator>(obsTerms[i]);
      for (auto x : pauliOp->getTerms().begin()->second.ops()) {
        if (x.second == "Z") {
          Z(x.first, i) = 1.0;
          nZeros++;
        } else if (x.second == "X") {
          X(x.first, i) = 1.0;
        } else if (x.second == "Y") {
          Z(x.first, i) = 1.0;
          X(x.first, i) = 1.0;
        }
      }

      if (nZeros == nOps)
        allZ.push_back(i);
    }

    // rankify if low-rank
    for (int i : allZ) {

      int j = 0;
      while (Z(j, i) != 1 && j < nOps) {
        j++;
      }

      Z.row(j).swap(X.row(j));
      auto H = provider->createInstruction("H", {(std::size_t)j});
      measurements->addInstruction(H);
    }

    // Perform Gauss-Jordan elimination on the X matrix
    // to bring it to the identity, while keeping track of the Z matrix
    //
    // Swapping rows
    for (int i = 0; i < nOps; i++) {

      if (X(i, i) == 0) {
        int c = 1;
        while ((i + c) < nOps && X(i + c, i) == 0)
          c++;
        if ((i + c) == nOps) {
          break;
        }
        X.row(i).swap(X.row(i + c));
        Z.row(i).swap(Z.row(i + c));

        // add a Swap gate if we swap rows
        auto swap = provider->createInstruction(
            "Swap", {(std::size_t)i, (std::size_t)(i + c)});
        measurements->addInstruction(swap);
      }

      for (int j = 0; j < nOps; j++) {

        if (i != j && X(j, i) != 0.0) {

          X.row(j) += X.row(i);
          Z.row(i) += Z.row(j);

          // enforcing modulo 2 addition
          for (int k = 0; k < nOps; k++) {
            if (X(j, k) == 2.0) {
              X(j, k) = 0.0;
            }
            if (Z(i, k) == 2.0) {
              Z(i, k) = 0.0;
            }
          }
          auto cx = provider->createInstruction(
              "CNOT", {(std::size_t)i, (std::size_t)j});
          measurements->addInstruction(cx);
        }
      }

      if (X == Eigen::MatrixXd::Identity(nOps, nOps))
        break;
    }

    // zeroing the 1's in the Z matrix
    for (int i = 0; i < nOps; i++) {

      // if 1 in the diagonal, add S
      if (Z(i, i) == 1.0) {
        auto s = provider->createInstruction("S", {(std::size_t)i});
        measurements->addInstruction(s);
        Z(i, i) = 0.0;
        continue;
      }

      // the Z matrix is symmetric, so if 1 in the off-diagonal,
      // add a CZ(i, j)
      for (int j = i + 1; j < nOps; j++) {
        if (Z(i, j) == 1.0) {
          auto cz = provider->createInstruction(
              "CZ", {(std::size_t)i, (std::size_t)j});
          measurements->addInstruction(cz);
          Z(i, j) = 0.0;
          Z(j, i) = 0.0;
          break;
        }
      }
    }

    // add H's and measure all qubits
    for (int i = 0; i < nOps; i++) {
      auto H = provider->createInstruction("H", {(std::size_t)i});
      measurements->addInstruction(H);
      auto measure = provider->createInstruction("Measure", {(std::size_t)i});
      measurements->addInstruction(measure);
    }
    groupMeasurements.push_back(measurements);
  }

  return groupMeasurements;
}

bool QubitWiseCommutativity::checkIfCommute(
    const std::map<int, std::string> op1,
    const std::map<int, std::string> op2) const {

  // first we check for common indices
  std::set<int> indices;
  for (auto &kv : op1) {
    for (auto &kvk : op2) {
      if (kv.first == kvk.first) {
        indices.insert(kv.first);
      }
    }
  }

  // if there are no common indices,
  // the operators commute because all
  // Paulis commmute with I
  if (indices.size() == 0) {
    return true;
  } else {

    // if there are common indices the operators qw commute
    // only if they are the same Pauli for the same index
    for (auto kv : indices) {
      if (op1.at(kv) != op2.at(kv)) {
        return false;
      }
    }
  }

  return true;
}

} // namespace xacc

REGISTER_PLUGIN(xacc::QubitWiseCommutativity, xacc::ObservableGrouping)