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
 *******************************************************************************/
#ifndef XACC_ALGORITHM_OPERATOR_POOLS_HPP_
#define XACC_ALGORITHM_OPERATOR_POOLS_HPP_

#include "adapt_vqe.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Observable.hpp"
#include "xacc_observable.hpp"
#include "FermionOperator.hpp"
#include "ObservableTransform.hpp"
//#include "<vector>"

namespace xacc{
namespace algorithm{

class UCCSD : public OperatorPool {

public:
  // check if pool exists
  bool isValidOperatorPool(const std::string &operatorPool) override {
    if(operatorPool == "UCCSD"){ // will add more pools in the future
      return true;
    }
    return false;
  }

  // generate the pool
  std::vector<std::shared_ptr<Observable>>
  generate(const int &nQubits, const int &nElectrons) override {

    const auto inv_sqrt_2 = "0.707106781186548";
    const auto one_inv_sqrt_12 = "0.288675134594813";
    const auto two_inv_sqrt_12 = "0.577350269189626";
    auto _nOccupied = (int)std::ceil(nElectrons / 2.0);
    auto _nVirtual = nQubits / 2 - _nOccupied;
    auto nSingle = _nOccupied * _nVirtual;
    auto nDouble = nSingle * (nSingle + 1) / 2;
    auto _nParameters = nSingle + nDouble;

    auto jw = xacc::getService<ObservableTransform>("jw");
    std::vector<std::shared_ptr<Observable>> pool;

    // single excitations
    for (int i; i < _nOccupied; i++){
      int ia = 2 * i;
      int ib = 2 * i + 1;
      for (int a; a < _nVirtual; a++){
        int aa = 2 * a + 2 * _nOccupied;
        int ab = 2 * a + 2 * _nOccupied + 1;

        // spin-adapted singles
        auto opString = constructOperatorString(aa, ia, inv_sqrt_2)
                      + constructOperatorString(ab, ib, inv_sqrt_2);
        auto fermiOp = xacc::quantum::getObservable("fermion", opString); 
        auto pauliOp = jw->transform(fermiOp);
        pool.push_back(pauliOp);
      }
    }

    // double excitations
    for (int i; i < _nOccupied; i++){
      int ia = 2 * i;
      int ib = 2 * i + 1;
      for (int j; j < _nOccupied; j++){
        int ja = 2 * j;
        int jb = 2 * j + 1;
        for (int a; a < _nVirtual; a++){
          int aa = 2 * a + 2 * _nOccupied;
          int ab = 2 * a + 2 * _nOccupied + 1;
          for ( int b; b < _nVirtual; b++){
            int ba = 2 * b + 2 * _nOccupied;
            int bb = 2 * b + 2 * _nOccupied + 1;

              auto opString = constructOperatorString(aa, ia, ba, ja, two_inv_sqrt_12)
                            + constructOperatorString(ab, ib, bb, jb, two_inv_sqrt_12)
                            + constructOperatorString(aa, ia, bb, jb, one_inv_sqrt_12)
                            + constructOperatorString(ab, ib, ba, ja, one_inv_sqrt_12)
                            + constructOperatorString(aa, ib, bb, ja, one_inv_sqrt_12)
                            + constructOperatorString(ab, ia, ba, jb, two_inv_sqrt_12);
              auto fermiOp = xacc::quantum::getObservable("fermion", opString);
              auto pauliOp = jw->transform(fermiOp);
              pool.push_back(pauliOp);

              opString = constructOperatorString(aa, ia, bb, jb, "0.5")
                        + constructOperatorString(ab, ib, ba, ja, "0.5")
                        + constructOperatorString(aa, ib, bb, ja, "0.5")
                        + constructOperatorString(ab, ia, ba, jb, "0.5");
              fermiOp = xacc::quantum::getObservable("fermion", opString);
              pauliOp = jw->transform(fermiOp);
              pool.push_back(pauliOp);
          }
        }
      }
    }
    return pool;
  }

  const std::string name() const override { return "uccsd"; }
  const std::string description() const override { return ""; }
};

} // namespace algorithm
} // namespace xacc

#endif
