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
#include <memory>
#include <string>

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
    
    auto _nOccupied = (int)std::ceil(nElectrons / 2.0);
    auto _nVirtual = nQubits / 2 - _nOccupied;
    auto _nOrbs = _nOccupied + _nVirtual;
    std::vector<std::shared_ptr<Observable>> pool;

    // single excitations
    for (int i = 0; i < _nOccupied; i++){
      int ia = i;
      int ib = i + _nOrbs;
      for (int a = 0; a < _nVirtual; a++){
        int aa = a + _nOccupied;
        int ab = a + _nOccupied + _nOrbs;

        // spin-adapted singles
        FermionOperator fermiOp;
        fermiOp = FermionOperator({{aa, 1}, {ia, 0}}, 0.5);
        fermiOp -= FermionOperator({{ia, 1}, {aa, 0}}, 0.5);

        fermiOp += FermionOperator({{ab, 1}, {ib, 0}}, 0.5);
        fermiOp -= FermionOperator({{ib, 1}, {ab, 0}}, 0.5);

        pool.push_back(std::dynamic_pointer_cast<Observable>(std::make_shared<FermionOperator>(fermiOp)));

      }
    }

    // double excitations
    for (int i = 0; i < _nOccupied; i++){
      int ia = i;
      int ib = i + _nOrbs;
      for (int j = i; j < _nOccupied; j++){
        int ja = j;
        int jb = j + _nOrbs;
        for (int a = 0; a < _nVirtual; a++){
          int aa = a + _nOccupied;
          int ab = a + _nOccupied + _nOrbs;
          for ( int b = a; b < _nVirtual; b++){
            int ba = b + _nOccupied;
            int bb = b + _nOccupied + _nOrbs;

            FermionOperator fermiOp;
            fermiOp = FermionOperator({{aa, 1}, {ia, 0}, {ba, 1}, {ja, 0}}, 2.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{ja, 1}, {ba, 0}, {ia, 1}, {aa, 0}}, 2.0 / std::sqrt(24.0));

            fermiOp += FermionOperator({{ab, 1}, {ib, 0}, {bb, 1}, {jb, 0}}, 2.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{jb, 1}, {bb, 0}, {ib, 1}, {ab, 0}}, 2.0 / std::sqrt(24.0));

            fermiOp += FermionOperator({{aa, 1}, {ia, 0}, {bb, 1}, {jb, 0}}, 1.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{jb, 1}, {bb, 0}, {ia, 1}, {aa, 0}}, 1.0 / std::sqrt(24.0));

            fermiOp += FermionOperator({{ab, 1}, {ib, 0}, {ba, 1}, {ja, 0}}, 1.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{ja, 1}, {ba, 0}, {ib, 1}, {ab, 0}}, 1.0 / std::sqrt(24.0));

            fermiOp += FermionOperator({{aa, 1}, {ib, 0}, {bb, 1}, {ja, 0}}, 1.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{ja, 1}, {bb, 0}, {ib, 1}, {aa, 0}}, 1.0 / std::sqrt(24.0));

            fermiOp += FermionOperator({{ab, 1}, {ia, 0}, {ba, 1}, {jb, 0}}, 1.0 / std::sqrt(24.0));
            fermiOp -= FermionOperator({{jb, 1}, {ba, 0}, {ia, 1}, {ab, 0}}, 1.0 / std::sqrt(24.0));

            pool.push_back(std::dynamic_pointer_cast<Observable>(std::make_shared<FermionOperator>(fermiOp)));

            fermiOp = FermionOperator({{aa, 1}, {ia, 0}, {bb, 1}, {jb, 0}}, 1.0 / (2.0 * std::sqrt(2.0)));
            fermiOp -= FermionOperator({{jb, 1}, {bb, 0}, {ia, 1}, {aa, 0}}, 1.0 / (2.0 * std::sqrt(2.0)));

            fermiOp += FermionOperator({{ab, 1}, {ib, 0}, {ba, 1}, {ja, 0}}, 1.0 / (2.0 * std::sqrt(2.0)));
            fermiOp -= FermionOperator({{ja, 1}, {ba, 0}, {ib, 1}, {ab, 0}}, 1.0 / (2.0 * std::sqrt(2.0)));

            fermiOp += FermionOperator({{aa, 1}, {ib, 0}, {bb, 1}, {ja, 0}}, -1.0 / (2.0 * std::sqrt(2.0)));
            fermiOp -= FermionOperator({{ja, 1}, {bb, 0}, {ib, 1}, {aa, 0}}, -1.0 / (2.0 * std::sqrt(2.0)));

            fermiOp += FermionOperator({{ab, 1}, {ia, 0}, {ba, 1}, {jb, 0}}, -1.0 / (2.0 * std::sqrt(2.0)));
            fermiOp -= FermionOperator({{jb, 1}, {ba, 0}, {ia, 1}, {ab, 0}}, -1.0 / (2.0 * std::sqrt(2.0)));

            pool.push_back(std::dynamic_pointer_cast<Observable>(std::make_shared<FermionOperator>(fermiOp)));

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