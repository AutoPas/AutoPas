/**
 * @file RegularGrid.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <memory>

#include "DomainDecomposition.h"
#include "autopas/AutoPas.h"
#include "mpi.h"
#include "src/TypeDefinitions.h"

class RegularGrid final : public DomainDecomposition {
 public:
  RegularGrid(int argc, char **argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
              const std::vector<double> &globalBoxMax);
  ~RegularGrid();

  using SharedAutoPasContainer = std::shared_ptr<autopas::AutoPas<ParticleType>>;

  void update() override;
  const int getDimensionCount() override { return _dimensionCount; }
  std::vector<double> getGlobalBoxMin() override { return _globalBoxMin; }
  std::vector<double> getGlobalBoxMax() override { return _globalBoxMax; }
  std::vector<double> getLocalBoxMin() override { return _localBoxMin; }
  std::vector<double> getLocalBoxMax() override { return _localBoxMax; }
  bool isInsideLocalDomain(const std::vector<double> &coordinates) override;

  int convertIdToIndex(const std::vector<int> &domainIndex);
  void exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer);
  void exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer);
  void receiveDataFromNeighbour(const int &neighbour, std::vector<char> &dataBuffer);
  void sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour);
  void synchronizeDomains();
  void waitForSendRequests();
  int getDomainIndex() { return _domainIndex; }

 private:
  // Global data
  int _dimensionCount;
  int _subdomainCount;
  std::vector<double> _globalBoxMin;
  std::vector<double> _globalBoxMax;
  std::vector<int> _decomposition;
  MPI_Comm _communicator;

  // Domain specific data
  int _domainIndex;
  std::vector<int> _domainId;
  std::vector<int> _neighbourDomainIndices;
  std::vector<double> _localBoxMin;
  std::vector<double> _localBoxMax;
  std::vector<MPI_Request> _sendRequests;
  std::vector<std::vector<char>> _sendBuffers;

  void initializeDecomposition();
  void initializeMPICommunicator();
  void initializeLocalDomain();
  void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
  void initializeLocalBox();
  void initializeNeighbourIds();

  void updateLocalBox();
  void sendParticles(std::vector<ParticleType> &particles, int &receiver);
  void receiveParticles(std::vector<ParticleType> &receivedParticles, int &source);
};
