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
#include "rbm_chemistry.hpp"

#include "AcceleratorBuffer.hpp"
#include "CompositeInstruction.hpp"
#include "Utils.hpp"
#include "xacc.hpp"
#include "xacc_service.hpp"
#include "xacc_observable.hpp"
#include "ObservableTransform.hpp"
#include "PauliOperator.hpp"
#include <cmath>
#include <map>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include "rbm.hpp"

#include <algorithm>
#include <vector>

using namespace xacc;
using namespace xacc::quantum;
// using namespace xacc::dwave;

using bitStringPair = std::pair<std::string, std::string>;
using Counts = std::map<std::string, int>;

namespace xacc {
namespace algorithm {
bool RBMChemistry::initialize(const HeterogeneousMap &parameters) {

  if (!parameters.pointerLikeExists<Observable>("observable")) {
    std::cout << "Obs was false\n";
    return false;
  }

  if (!parameters.keyExists<int>("hidden-layers")) {
    std::cout << "Number of hidden layers was false\n";
    return false;
  }

  if (!parameters.pointerLikeExists<Accelerator>("accelerator")) {
    std::cout << "Acc was false\n";
    return false;
  }

  observable = parameters.getPointerLike<Observable>("observable");
  accelerator = parameters.getPointerLike<Accelerator>("accelerator");
  nHidden = parameters.get<int>("hidden-layers");
  nVisible = observable->nBits();

  // use predefined RBM
  rbm = std::dynamic_pointer_cast<CompositeInstruction>(
      xacc::getService<Instruction>("rbm-ap"));
  rbm->expand({{"nv", nVisible}, {"nh", nHidden}});

  // Get sparse matrix for Hamiltonian
  auto JWObservable = std::dynamic_pointer_cast<PauliOperator>(
      xacc::getService<ObservableTransform>("jw")->transform(
          xacc::as_shared_ptr(observable)));
  auto sparseObservable = JWObservable->to_sparse_matrix();
  for (auto ss : sparseObservable) {
    auto bra = generateBitString(ss.row());
    auto ket = generateBitString(ss.col());
    auto statePair = std::make_pair(bra, ket);
    hamiltonian.push_back(std::make_pair(statePair, std::real(ss.coeff())));
  }

  return true;
}

const std::vector<std::string> RBMChemistry::requiredParameters() const {
  return {"observable", "hidden-layers", "accelerator"};
}

void RBMChemistry::execute(
    const std::shared_ptr<AcceleratorBuffer> buffer) const {

  // store connections
  std::vector<std::pair<int, int>> connections;
  for (int i = 0; i < nVisible; i++) {
    for (int j = 0; j < nHidden; j++) {
      connections.push_back(std::make_pair(i, j + nVisible));
    }
  }

  // randomly initialize weights and biases
  // std::vector<double> weights_biases =
  //    random_vector(-wInitMax, wInitMax,
  //                  nVisible + nHidden + 1 + nVisible * nHidden + nVisible);

  // this is the set of pseudo random numbers from the seed used in
  // https://code.ornl.gov/aqw/rbm-chemistry/-/blob/master/sindhu/methods/gd.py
  std::vector<double> weights_biases = {
      0.09340596780273533,   0.009446449835144463, 0.09453687199297686,
      0.04296319873487292,   0.03954576491945416,  -0.05678210088392473,
      0.09525489095524836,   -0.09875394895908203, -0.04940352752331121,
      -0.013041693519110845, 0.05587658435875051,  -0.060462985079949384,
      0.07259864711984446,   0.09668013543506257,  -0.06723155171906026,
      0.019466788786571848,  -0.098202780466489,   -0.022685743471274125,
      -0.0911679884137001,   0.0913305935428472,   -0.01277067062404047,
      0.08979546135631256,   0.0572611971870122,   0.07325785971633969,
      -0.06536691570081038,  -0.0850102825973786,  0.020148544275553862,
      -0.06640556325631924,  0.04667603350211397};

  // pass weights and biases to RBM
  auto program = rbm->operator()(std::vector<double>(
      weights_biases.begin(),
      weights_biases.begin() + nVisible + nHidden + nVisible * nHidden));

  // Sample visible states
  auto tmpBuffer = qalloc(buffer->size());
  auto counts = sample(tmpBuffer, program);

  // Update parameters
  Eigen::VectorXd update = computeUpdate(weights_biases, counts);
}

// This follows
// https://code.ornl.gov/aqw/rbm-chemistry/-/blob/master/sindhu/methods/gd.py#L54
Eigen::VectorXd RBMChemistry::computeUpdate(std::vector<double> wb,
                                            Counts counts) const {

  // assign weights and biases
  auto a = std::vector<double>(wb.begin(), wb.begin() + nVisible);
  auto b = std::vector<double>(wb.begin() + nVisible,
                               wb.begin() + nVisible + nHidden);
  auto w =
      std::vector<double>(wb.begin() + nVisible + nHidden,
                          wb.begin() + nVisible + nHidden + nVisible * nHidden);
  auto d = std::vector<double>(
      wb.begin() + nVisible + nHidden + nVisible * nHidden, wb.end() - 1);
  auto c = wb.back();

  std::vector<std::string> bt = {
      "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
      "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111",
  };

  Counts visibleCounts;
  auto nStates = (int)std::pow(2, nVisible);
  Eigen::VectorXd totalCounts = Eigen::VectorXd::Zero(nStates);
  for (auto k : counts) {
    visibleCounts.emplace(k.first.substr(0, nVisible), k.second);
  }

  // add all 2^nVisible states with 0 counts
  for (int i = 0; i < nStates; i++) {

    if (visibleCounts.find(generateBitString(i)) == visibleCounts.end()) {
      visibleCounts.emplace(generateBitString(i), 0);
    } else {
      totalCounts(i) = (double)visibleCounts[generateBitString(i)];
      continue;
    }
  }

  Eigen::MatrixXd states(nStates, nVisible);
  Eigen::VectorXd d_ = Eigen::VectorXd::Map(d.data(), d.size());
  Eigen::VectorXd b_ = Eigen::VectorXd::Map(b.data(), b.size());

  // Get bit strings and map bits to 1/-1
  for (int i = 0; i < nStates; i++) {
    auto bitString = generateBitString(i);
    for (int j = 0; j < nVisible; j++) {
      if (bitString[j]) {
        states(i, j) = 1;
      } else {
        states(i, j) = -1;
      }
    }
  }

  Eigen::VectorXd stateCounts = totalCounts / totalCounts.sum();

  // Signs
  Eigen::VectorXd sx(nStates);
  for (int i = 0; i < nStates; i++) {
    sx(i) = tanh(states.row(i) * d_ + c);
  }

  // Compute probabilities
  Eigen::VectorXd probabilities(nStates);
  for (int i = 0; i < nStates; i++) {
    probabilities(i) = (std::pow(sx(i), 2) * std::pow(stateCounts(i), 3));
  }

  auto totalProbability = probabilities.sum();
  for (int i = 0; i < nStates; i++) {
    probabilities(i) = probabilities(i) / totalProbability;
  }

  // Psi
  Eigen::VectorXd psi(nStates);
  for (int i = 0; i < nStates; i++) {
    psi(i) = sqrt(probabilities(i));
    if (sx(i) < 0.0)
      psi(i) *= -1.0;
  }

  // Energies
  Eigen::VectorXd localEnergies(nStates);
  for (int i = 0; i < nStates; i++) {
    double localEnergy = 0.0;
    for (int j = 0; j < nStates; j++) {

      bitStringPair braket(generateBitString(i), generateBitString(j));
      for (auto term : hamiltonian)
        if (braket == term.first) {
          localEnergy += term.second * psi(j);
          break;
        }
    }
    localEnergies(i) = (localEnergy / psi(i));
  }

  double energy = 0.0;
  for (int i = 0; i < nStates; i++) {
    energy += localEnergies(i) * probabilities(i);
  }

  // I may come back here to add the best set of parameters
  Eigen::MatrixXd dA = 3 * states.transpose() / 2.0;
  Eigen::Map<Eigen::MatrixXd> w_(w.data(), nVisible, nHidden);
  Eigen::MatrixXd bMatrix(nHidden, nStates);
  for (int i = 0; i < nStates; i++)
    bMatrix.col(i) = b_;
  Eigen::MatrixXd theta = w_ * states.transpose() + bMatrix;
  Eigen::MatrixXd dB =
      3 * theta.unaryExpr([](double x) { return tanh(x) / 2.0; });

  Eigen::MatrixXd dW(nVisible * nHidden, nStates);
  for (int i = 0; i < nVisible; i++) {
    for (int j = 0; j < nHidden; j++) {
      dW.row(i * nHidden + j) =
          states.transpose().row(i).cwiseProduct(dB.row(j));
    }
  }

  Eigen::VectorXd sxInverse(sx.size());
  for (int i = 0; i < sx.size(); i++) {
    auto sInverse = 1.0 / sx(i);
    if (sInverse <= -1.0e11) {
      sxInverse(i) = -1.0e11;
    } else if (sInverse >= 1.0e11) {
      sxInverse(i) = 1.0e11;
    } else {
      sxInverse(i) = sInverse;
    }
  }
  Eigen::VectorXd dc = (sxInverse - sx); //.transpose();

  Eigen::MatrixXd dd(nVisible, nStates);
  for (int i = 0; i < nVisible; i++) {
    dd.row(i) = states.transpose().row(i).cwiseProduct(dc.transpose());
  }

  Eigen::MatrixXd dP(dA.rows() + dB.rows() + dW.rows() + dd.rows() +
                         dc.transpose().rows(),
                     dA.cols());
  dP << dA, dB, dW, dd, dc.transpose();

  Eigen::MatrixXd E_dP = dP * probabilities;
  Eigen::MatrixXd dP_conj = dP.conjugate();
  Eigen::MatrixXd E_dP_conj = dP_conj * probabilities;
  Eigen::MatrixXd dP_probs(dP.rows(), dP.cols());
  for (int i = 0; i < dP.rows(); i++) {
    dP_probs.row(i) = dP.row(i).cwiseProduct(probabilities.transpose());
  }

  Eigen::MatrixXd S =
      dP_conj * dP_probs.transpose() - E_dP_conj * E_dP.transpose();

  Eigen::MatrixXd F =
      dP_conj * (localEnergies.array() * probabilities.array()).matrix() -
      energy * E_dP_conj;

  Eigen::VectorXd update =
      (S + eps * Eigen::MatrixXd::Identity(S.rows(), S.cols()))
          .colPivHouseholderQr()
          .solve(F);

  return update;
}

std::vector<double> RBMChemistry::random_vector(const double l_range,
                                                const double r_range,
                                                const std::size_t size) const {
  // Generate a random initial parameter set
  std::random_device rnd_device;
  std::mt19937 mersenne_engine{rnd_device()}; // Generates random integers
  std::uniform_real_distribution<double> dist{l_range, r_range};
  auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
  std::vector<double> vec(size);
  std::generate(vec.begin(), vec.end(), gen);
  return vec;
}

// generates a bit string from an integer
std::string RBMChemistry::generateBitString(unsigned decimal) const {
  std::string bitString;
  do {
    bitString.push_back('0' + (decimal & 1));
  } while (decimal >>= 1);
  std::reverse(bitString.begin(), bitString.end());
  return bitString;
}

// This is for run on actual hardware
double RBMChemistry::computeScale(
    const std::vector<double> a, const std::vector<double> b,
    const std::vector<double> w,
    const std::map<int, std::vector<int>> embedding) const {

  std::vector<double> a_(a.size()), b_(a.size());
  for (int i = 0; i < a.size(); i++) {
    if (embedding.find(i) != embedding.end()) {
      a_[i] = fabs(a[i] / embedding.find(i)->second.size());
    }
  }

  for (int i = 0; i < b.size(); i++) {
    if (embedding.find(i) != embedding.end()) {
      b_[i] = fabs(b[i] / embedding.find(i)->second.size());
    }
  }

  auto upperLimit =
      std::min(2.0 / *std::max_element(a_.begin(), a_.end()),
               std::min(*std::max_element(b_.begin(), b_.end()),
                        1.0 / *std::max_element(w.begin(), w.end())));

  auto lowerLimit = minCoefficient /
                    std::min(*std::min_element(a_.begin(), a_.end()),
                             std::min(*std::min_element(b_.begin(), b_.end()),
                                      *std::min_element(w.begin(), w.end())));

  if (lowerLimit > upperLimit) {
    xacc::warning(
        "Lower limit greater than upper limit while computing scale.");
    return upperLimit;
  } else {
    return lowerLimit;
  }
}

/*
void RBMChemistry::getEmbedding(
    const std::shared_ptr<CompositeInstruction> problem,
    std::shared_ptr<AcceleratorBuffer> &buffer) const {

  std::map<int, std::vector<int>> embedding;
  auto hardwareEdges = accelerator->getConnectivity();
  auto probGraph = problem->toGraph();
  int maxBit = 0;
  for (auto &e : hardwareEdges) {
    if (e.first > maxBit) {
      maxBit = e.first;
    }
    if (e.second > maxBit) {
      maxBit = e.second;
    }
  }

  auto hardwareGraph = xacc::getService<Graph>("boost-ugraph");
  for (int i = 0; i < maxBit + 1; i++) {
    HeterogeneousMap props{std::make_pair("bias", 1.0)};
    hardwareGraph->addVertex(props); //
  }

  for (auto &e : hardwareEdges) {
    hardwareGraph->addEdge(e.first, e.second);
  }

  if (!buffer->hasExtraInfoKey("embedding")) {

    auto embeddingAlgo = xacc::getService<EmbeddingAlgorithm>(default_emb_algo);
    embedding = embeddingAlgo->embed(probGraph, hardwareGraph);
    buffer->addExtraInfo("embedding", embedding);

  } else {
    auto tmppp = buffer->getInformation("embedding")
                     .as<std::map<int, std::vector<int>>>();
    for (auto &kv : tmppp) {
      embedding.insert(kv);
    }
  }
}
*/

Counts
RBMChemistry::sample(std::shared_ptr<AcceleratorBuffer> &buffer,
                     std::shared_ptr<CompositeInstruction> program) const {

  if (accelerator->name() == "dwave-neal") {

    // neal simulator is the default
    auto q = qalloc(nVisible + nHidden);
    accelerator->execute(q, program);
    return q->getMeasurementCounts();

  } else {

    // this is for actual DWave hardware

    // first determine the embedding
    std::map<int, std::vector<int>> embedding;
    auto hardwareEdges = accelerator->getConnectivity();
    auto probGraph = program->toGraph();
    int maxBit = 0;
    for (auto &e : hardwareEdges) {
      if (e.first > maxBit) {
        maxBit = e.first;
      }
      if (e.second > maxBit) {
        maxBit = e.second;
      }
    }

    auto hardwareGraph = xacc::getService<Graph>("boost-ugraph");
    for (int i = 0; i < maxBit + 1; i++) {
      HeterogeneousMap props{std::make_pair("bias", 1.0)};
      hardwareGraph->addVertex(props); //
    }

    for (auto &e : hardwareEdges) {
      hardwareGraph->addEdge(e.first, e.second);
    }

    if (!buffer->hasExtraInfoKey("embedding")) {

      auto embeddingAlgo =
          xacc::getService<EmbeddingAlgorithm>(default_emb_algo);
      embedding = embeddingAlgo->embed(probGraph, hardwareGraph);
      buffer->addExtraInfo("embedding", embedding);

    } else {

      auto tmppp = buffer->getInformation("embedding")
                       .as<std::map<int, std::vector<int>>>();
      for (auto &kv : tmppp) {
        embedding.insert(kv);
      }
    }

    return {};
  }

  // TODO: need to finish this
  // and add an option to simulate using other backends (TNQVM/ExaTN)
}

} // namespace algorithm
} // namespace xacc
