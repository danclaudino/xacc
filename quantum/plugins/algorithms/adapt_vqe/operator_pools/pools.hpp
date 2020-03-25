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

  UCCSD() = default;

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
    
    const double inv_sqrt_2 = 0.707106781186548;
    const double one_inv_sqrt_12 = 0.288675134594813;
    const double two_inv_sqrt_12 = 0.577350269189626;
    auto _nOccupied = (int)std::ceil(nElectrons / 2.0);
    auto _nVirtual = nQubits / 2 - _nOccupied;

    auto jw = xacc::getService<ObservableTransform>("jw");
    std::vector<std::shared_ptr<Observable>> pool;

    // single excitations
    for (int i = 0; i < _nOccupied; i++){
      int ia = 2 * i;
      int ib = 2 * i + 1;
      for (int a = 0; a < _nVirtual; a++){
        int aa = 2 * a + 2 * _nOccupied;
        int ab = 2 * a + 2 * _nOccupied + 1;

        // spin-adapted singles
        FermionOperator antiHermitianOp;
        antiHermitianOp = FermionOperator({{aa, 1}, {ia, 0}}, inv_sqrt_2);
        antiHermitianOp -= FermionOperator({{ia, 1}, {aa, 0}}, inv_sqrt_2);

        antiHermitianOp += FermionOperator({{ab, 1}, {ib, 0}}, inv_sqrt_2);
        antiHermitianOp -= FermionOperator({{ib, 1}, {ab, 0}}, inv_sqrt_2);

        // std::cout << antiHermitianOp.toString();
        auto pauliOp = jw->transform(std::shared_ptr<Observable>(&antiHermitianOp, [](Observable *) {}));
        pool.push_back(pauliOp);
      }
    }

    // double excitations
    for (int i = 0; i < _nOccupied; i++){
      int ia = 2 * i;
      int ib = 2 * i + 1;
      for (int j = i; j < _nOccupied; j++){
        int ja = 2 * j;
        int jb = 2 * j + 1;
        for (int a = 0; a < _nVirtual; a++){
          int aa = 2 * a + 2 * _nOccupied;
          int ab = 2 * a + 2 * _nOccupied + 1;
          for ( int b = a; b < _nVirtual; b++){
            int ba = 2 * b + 2 * _nOccupied;
            int bb = 2 * b + 2 * _nOccupied + 1;

            FermionOperator antiHermitianOp;
            antiHermitianOp = FermionOperator({{aa, 1}, {ia, 0}, {ba, 1}, {ja, 0}}, two_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{ja, 1}, {ba, 0}, {ia, 1}, {aa, 0}}, two_inv_sqrt_12);

            antiHermitianOp += FermionOperator({{ab, 1}, {ib, 0}, {bb, 1}, {jb, 0}}, two_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{jb, 1}, {bb, 0}, {ib, 1}, {ab, 0}}, two_inv_sqrt_12);

            antiHermitianOp += FermionOperator({{aa, 1}, {ia, 0}, {bb, 1}, {jb, 0}}, one_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{jb, 1}, {bb, 0}, {ia, 1}, {aa, 0}}, one_inv_sqrt_12);

            antiHermitianOp += FermionOperator({{ab, 1}, {ib, 0}, {ba, 1}, {ja, 0}}, one_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{ja, 1}, {ba, 0}, {ib, 1}, {ab, 0}}, one_inv_sqrt_12);

            antiHermitianOp += FermionOperator({{aa, 1}, {ib, 0}, {bb, 1}, {ja, 0}}, one_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{ja, 1}, {bb, 0}, {ib, 1}, {aa, 0}}, one_inv_sqrt_12);

            antiHermitianOp += FermionOperator({{ab, 1}, {ia, 0}, {ba, 1}, {jb, 0}}, one_inv_sqrt_12);
            antiHermitianOp -= FermionOperator({{jb, 1}, {ba, 0}, {ia, 1}, {ab, 0}}, one_inv_sqrt_12);

            auto pauliOp = jw->transform(std::shared_ptr<Observable>(&antiHermitianOp, [](Observable *) {}));
            pool.push_back(pauliOp);

            antiHermitianOp = FermionOperator({{aa, 1}, {ia, 0}, {bb, 1}, {jb, 0}}, 0.5);
            antiHermitianOp -= FermionOperator({{jb, 1}, {bb, 0}, {ia, 1}, {aa, 0}}, 0.5);

            antiHermitianOp += FermionOperator({{ab, 1}, {ib, 0}, {bb, 1}, {jb, 0}}, 0.5);
            antiHermitianOp -= FermionOperator({{jb, 1}, {bb, 0}, {ib, 1}, {ab, 0}}, 0.5);

            antiHermitianOp += FermionOperator({{aa, 1}, {ib, 0}, {bb, 1}, {ja, 0}}, -0.5);
            antiHermitianOp -= FermionOperator({{ja, 1}, {bb, 0}, {ib, 1}, {aa, 0}}, -0.5);

            antiHermitianOp = FermionOperator({{ab, 1}, {ia, 0}, {ba, 1}, {jb, 0}}, 0.5);
            antiHermitianOp -= FermionOperator({{jb, 1}, {ba, 0}, {ia, 1}, {ab, 0}}, 0.5);

            pauliOp = jw->transform(std::shared_ptr<Observable>(&antiHermitianOp, [](Observable *) {}));
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
