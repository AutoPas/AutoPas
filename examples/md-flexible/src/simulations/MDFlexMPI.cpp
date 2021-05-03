/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"


// @todo remove as this is currently only used for debugging
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ArrayMath.h"

#include <algorithm>
#include <iostream>

namespace {
	void sendParticles(std::vector<Particle> &particles, int receiver, MPI_Comm &communicator, MPI_Request &sendRequest, std::vector<double> &buffer) {
  	for (auto &particle : particles) {
    	std::vector<double> serializedParticle = particle.serialize();
    	buffer.insert(std::end(buffer), std::begin(serializedParticle), std::end(serialized));
  	}
  	MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, 0, comm, &sendRequest);
	}

	void receiveParticles(std::vector<Particle> &receivedParticles, int sender, MPI_Comm &communicator, MPI_Request &receiveRequest) {
		// @todo: finish implementation of this function
  	MPI_Status status;
  	MPI_Probe(sender, 0, communicator, &status);

  	int receiveBufferSize;
  	MPI_Get_count(&status, MPI_DOUBLE, &receiveBufferSize);
  	std::vector<double> receiveBuffer(receiveBufferSize);

  	MPI_Recv(receiveBuffer.data(), receiveBufferSize, MPI_DOUBLE, neighbor, 0, communicator, MPI_STATUS_IGNORE);

		// @todo: serialize and deserialize are not defined. Either create an new MPI_Datatyep or implement these functions for the particle
  	for (size_t i = 0; i < (size_t)receiveBufferSize;) {
    	auto particle = Particle::deserialize(receiveBuffer.data(), i);
    	receiveParticles.push_back(particle);
  	}
	}
}

MDFlexMPI::MDFlexMPI(int argc, char** argv) : MDFlexSimulation(argc, argv){
  MPI_Init(&argc, &argv);

	int dimensionCount = _configuration->boxMin.value.size();
	double* globalBoxMin = &(_configuration->boxMin.value[0]);
	double* globalBoxMax = &(_configuration->boxMax.value[0]);
	_domainDecomposition = std::make_unique<RegularGridDecomposition>(dimensionCount, globalBoxMin, globalBoxMax);

  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(_configuration->cutoff.value);

	std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();

	for (int i = 0; i < dimensionCount; ++i){
		_configuration->boxMin.value[i] = localBoxMin[i];
		_configuration->boxMax.value[i] = localBoxMax[i];
	}

	_autoPasContainer = std::make_unique<autopas::AutoPas<Simulation::ParticleType>>(std::cout);	
}

MDFlexMPI::~MDFlexMPI() {
	MPI_Finalize();
}

void MDFlexMPI::run(){

	// @todo: make variable part of MDFlexConfig
	int iterationsPerSuperstep = 10;
	int remainingIterations = _configuration->iterations.value;

	for (int i = 0; i < _configuration->iterations.value; i+=iterationsPerSuperstep){
	//	_simulation->initialize(*_configuration, autopasContainer);
	//_simulation->simulate(autopasContainer);

		executeSuperstep(iterationsPerSuperstep);

	//	remainingIterations -= iterationsPerSuperstep;
	}

	//_simulation->printStatistics(autopasContainer);
}

void MDFlexMPI::executeSuperstep(const int iterationsPerSuperstep){
	// Update particles
	updateParticles(iterationsPerSuperstep);

	MPI_Barrier(_domainDecomposition->getCommunicator());

	// Update AutoPas container
	auto [emigrants, updated] = _autoPasContainer->updateContainer();
	if (updated) {
 		// @todo: implement this function using the 2-step scheme
		sendEmigrantsToNeighbours(emigrants);
	}
}

void MDFlexMPI::updateParticles(const int iterationsPerSuperstep){
	double deltaT = _configuration->deltaT.value;
	for (int i = 0; i < iterationsPerSuperstep; ++i){
 		for (auto particle = _autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
   		auto v = particle->getV();
   		auto m = _particlePropertiesLibrary->getMass(particle->getTypeId());
   		auto f = particle->getF();
   		particle->setOldF(f);
   		particle->setF({0., 0., 0.});
   		v = autopas::utils::ArrayMath::mulScalar(v, deltaT);
   		f = autopas::utils::ArrayMath::mulScalar(f, (deltaT * deltaT / (2 * m)));
   		auto newR = autopas::utils::ArrayMath::add(v, f);
   		particle->addR(newR);
 		}
	}
}

void MDFlexMPI::sendEmigrantsToNeighbours(std::vector<Particle> &emigrants){
	std::vector<double> localBoxMin = _decomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _decomposition->getLocalBoxMax();
	int dimensionCount = localBoxMin.size();

	MPI_Comm communicator = _decomposition->getCommunicator();

	std::vector<Particle> immigrant;
	std::vector<Particle> receivedParticles = emigrants;

	std::vector<Particle> particlesForPreceedingNeighbour;
	std::vector<Particle> particlesForSucceedingNeighbour;

	std::array<MPI_Request, 4> mpiRequests;
	std::array<std::vector<double>, 2> sendBuffers;

	for(int i = 0; i < dimensionCount; ++i){
		for (auto &particle : emigrans){
			if (particle->GetR()[i] - localBoxMin[i] < 0){
				particlesForPreceedingNeighbour.push_back(particle.get());
			}
			else (particle->GetR()[i] - localBoxMax[i] > 0 {
				particlesForSucceedingNeighbour.push_back(particle.get());
			}
		}

		sendParticles(particlesForPreceedingNeighbour.data(), dimensionCount * 2, communicator, &mpiRequests[0], &sendBuffers[0]);
		sendParticles(particlesForPreceedingNeighbour.data(), dimensionCount * 2 + 1, communicator, &mpiRequests[1], &sendBuffers[1]);

		// @todo: Receive particles from preceeding neighbour
		// @todo: Receive particles from succeeding neighbour

		// @todo: replace with MPI_wait
		MPI_Barrier(communicator);

		particlesForPreceedingNeighbour.clear();
		particlesForSucceedingNeighbour.clear();
		sendBuffer[0].clear();
		sendBuffer[1].clear();

		int nextDimensionIndex = (i + 1) % dimensionCount;
		for (auto &particle : receivedParticles){
			if (particle->GetR()[nextDimensionIndex] - localBoxMin[nextDimensionIndex] < 0){
				particlesForPreceedingNeighbour.push_back(particle.get());
			}
			else if (particle->GetR()[nextDimensionIndex] - localBoxMax[nextDimensionIndex] > 0 {
				particlesForSucceedingNeighbour.push_back(particle.get());
			}
			else {
				immigratingParticles.push_back(particle.get());
			}
		}

		// @todo: send migrating paricles to preceeding neighbour in next dimension
		// @todo: send migrating particles to succeeding neighbour in next dimension
		// @todo: receive migrating particles from preceeding neighbour in next dimension
		// @todo: receive migrating particles from succeeding neighbour in next dimension

		MPI_Barrier(communicator);

		particlesForPreceedingNeighbour.clear();
		particlesForSucceedingNeighbour.clear();
	}

	for (auto &particle : imigratingParticles) {
		_autoPasContainer->addParticle(particle);
	}
}
