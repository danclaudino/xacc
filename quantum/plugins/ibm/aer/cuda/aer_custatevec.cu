#include "aer_custatevec.hpp"
#include "controllers/aer_controller.hpp"

using StateVector = AER::Statevector::State<AER::QV::QubitVectorThrust<double>>;

void* getState() {
  StateVector* state = new StateVector();
  return state;
}
