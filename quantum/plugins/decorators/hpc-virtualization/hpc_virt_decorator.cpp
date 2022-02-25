/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *   Daniel Claudino - MPI native implementation
 *******************************************************************************/
#include "hpc_virt_decorator.hpp"
#include "InstructionIterator.hpp"
#include "Utils.hpp"
#include "xacc.hpp"
#include <numeric>

namespace xacc {
namespace quantum {

void HPCVirtDecorator::initialize(const HeterogeneousMap &params) {
  decoratedAccelerator->initialize(params);

  if (params.keyExists<int>("n-virtual-qpus")) {
    if (qpuComm && n_virtual_qpus != params.get<int>("n-virtual-qpus")) {
      // We don't support changing the number of virtual QPU's
      // i.e. between xacc::Initialize and xacc::Finalize,
      // we must use an HPCVirtDecorator with a consistent number of virtual
      // QPU's.
      xacc::error(
          "Dynamically changing the number of virtual QPU's is not supported.");
    }
    n_virtual_qpus = params.get<int>("n-virtual-qpus");
  }
}

void HPCVirtDecorator::updateConfiguration(const HeterogeneousMap &config) {
  decoratedAccelerator->updateConfiguration(config);
}

void HPCVirtDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> function) {

  if (decoratedAccelerator)
    decoratedAccelerator->execute(buffer, function);

  return;
}

void HPCVirtDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>> functions) {

  // The goal here is the following:
  // Assume we have an HPC system with a set of
  // compute ranks M, and we partition that set of ranks into
  // communicator sub-groups, each with a dedicated virtual QPU available.
  // Also we are given a vector of observable terms
  // (vec of CompositeInstructions) of size N (N>>1),
  // we wish to evaluate the <O> for each by partitioning
  // their quantum execution across the node sub-groups.


  // Get the rank and size in the original communicator
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  auto world = std::make_shared<ProcessGroup>(MPI_COMM_WORLD, world_size);

  if (world_size < n_virtual_qpus) {
    // The number of MPI processes is less than the number of requested virtual
    // QPUs, just execute as if there is only one virtual QPU and give the QPU
    // the whole MPI_COMM_WORLD.
    void *qpu_comm_ptr = reinterpret_cast<void *>(world.get());
    if (!qpuComm) {
      qpuComm = world;
    }
    decoratedAccelerator->updateConfiguration(
        {{"mpi-communicator", qpu_comm_ptr}});
    // just execute
    decoratedAccelerator->execute(buffer, functions);
    return;
  }

  // get the number of sub-communicators
  // Everybody split the CompositeInstructions vector into n_virtual_qpu segments
  // and get the color for the given process
  int color;
  auto split_vecs = split_vector(functions, n_virtual_qpus, world_rank, color);

  // Split the communicator based on the color and use the
  // original rank for ordering
  if (!qpuComm) {
    // Splits MPI_COMM_WORLD into sub-communicators if not already.
    qpuComm = world->split(color);
  }

  // current rank now has a color to indicate which sub-comm it belongs to
  // Give that sub communicator to the accelerator
  void *qpu_comm_ptr = reinterpret_cast<void *>(qpuComm->getMPICommProxy().getRef<MPI_Comm>());
  decoratedAccelerator->updateConfiguration(
      {{"mpi-communicator", qpu_comm_ptr}});

  // Get the segment corresponding to this color
  auto my_circuits = split_vecs[color];

  // Create a local buffer and execute
  auto my_buffer = xacc::qalloc(buffer->size());
  decoratedAccelerator->execute(my_buffer, my_circuits);

  // Create a mapping of all local buffer names to the buffer
  std::map<std::string, std::shared_ptr<AcceleratorBuffer>>
      local_name_to_buffer;
  for (auto &child : my_buffer->getChildren()) {
    local_name_to_buffer.insert({child->name(), child});
  }

  // Every sub-group compute local energy
  double local_energy = 0.0;
  for (auto &my_circuit : my_circuits) {
    local_energy +=
        std::real(my_circuit->getCoefficient()) *
        local_name_to_buffer[my_circuit->name()]->getExpectationValueZ();
  }

  // Split world along rank-0 in each sub-communicator and reduce the local energies
  auto zeroRanksComm = world->split(world_rank == qpuComm->getProcessRanks()[0]);
  double global_energy = 0.0;
  if (world_rank == qpuComm->getProcessRanks()[0]) {
    MPI_Allreduce(&local_energy, &global_energy, 1,
      MPI_DOUBLE, MPI_SUM, zeroRanksComm->getMPICommProxy().getRef<MPI_Comm>());
  }

  MPI_Barrier(qpuComm->getMPICommProxy().getRef<MPI_Comm>());  

  // Now broadcast that global_energy to all
  // other ranks in the sub groups
  MPI_Bcast(&global_energy, 1, MPI_DOUBLE, 0, qpuComm->getMPICommProxy().getRef<MPI_Comm>());

  // Setup a barrier
  MPI_Barrier(qpuComm->getMPICommProxy().getRef<MPI_Comm>());  
  MPI_Barrier(world->getMPICommProxy().getRef<MPI_Comm>());  

  // every rank should have global_energy now
  buffer->addExtraInfo("__internal__decorator_aggregate_vqe__", global_energy);
  buffer->addExtraInfo("rank", world_rank);

  return;
}

} // namespace quantum
} // namespace xacc

#include "cppmicroservices/BundleActivator.h"
#include "cppmicroservices/BundleContext.h"
#include "cppmicroservices/ServiceProperties.h"

using namespace cppmicroservices;

namespace {

class US_ABI_LOCAL HPCVirtActivator : public BundleActivator {

public:
  HPCVirtActivator() {}

  void Start(BundleContext context) {
    auto c = std::make_shared<xacc::quantum::HPCVirtDecorator>();

    context.RegisterService<xacc::AcceleratorDecorator>(c);
    context.RegisterService<xacc::Accelerator>(c);

  }

  /**
   */
  void Stop(BundleContext /*context*/) {}
};

} // namespace

CPPMICROSERVICES_EXPORT_BUNDLE_ACTIVATOR(HPCVirtActivator)
