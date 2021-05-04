/**
 * @file MDFlexMPI.cpp
 * @author J. Körner
 * @date 07.04.2021
 */

#include "MDFlexMPI.h"

#include "autopas/utils/ArrayMath.h"

#include <algorithm>
#include <iostream>

namespace {
	struct Attributes {
		// ParticleBase attributes
		double positionX; 		
		double positionY;
		double positionZ;
		double velocityX; 		
		double velocityY;
		double velocityZ;
		double forceX;
		double forceY;
		double forceZ;
 		unsigned long id; 		

		// MoleculeLJ attributes
  	size_t typeId;
		double oldForceX;
		double oldForceY;
		double oldForceZ;
	};

	void transferParticleAttributes(ParticleType &particle, Attributes &attributes){
		attributes.positionX = particle.getR()[0];
		attributes.positionY = particle.getR()[1];
		attributes.positionZ = particle.getR()[2];

		attributes.velocityX = particle.getV()[0];
		attributes.velocityY = particle.getV()[1];
		attributes.velocityZ = particle.getV()[2];

		attributes.forceX = particle.getF()[0];
		attributes.forceY = particle.getF()[1];
		attributes.forceZ = particle.getF()[2];

		attributes.id = particle.getID();

		attributes.typeId = particle.getTypeId();
		attributes.oldForceX = particle.getOldf()[0];
		attributes.oldForceY = particle.getOldf()[1];
		attributes.oldForceZ = particle.getOldf()[2];
	}

	void applyAttributestoParticle(Attributes &attributes, ParticleType &particle){
		// @todo: implement this function
	}

	const char* serializeAttributes(ParticleType &particle, Attributes &attributes){
 		return reinterpret_cast<const char*>(&(particle.getR()[0])); 
	}

	void sendParticles(std::vector<ParticleType> &particles, int &receiver, MPI_Comm &communicator, MPI_Request &sendRequest, std::vector<char> &buffer) {
  	for (auto &particle : particles) {
			Attributes attributes;
			transferParticleAttributes(particle, attributes);

			const char* attributesStart = serializeAttributes(particle, attributes);
			const char* attributesEnd = attributesStart + sizeof(Attributes) - 1;

    	buffer.insert(std::end(buffer), attributesStart, attributesEnd);
  	}

  	MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, receiver, 0, communicator, &sendRequest);
	}

	void receiveParticles(std::vector<ParticleType> &receivedParticles, int source, MPI_Comm &communicator) {
  	MPI_Status status;
  	MPI_Probe(source, 0, communicator, &status);

  	int receiveBufferSize;
  	MPI_Get_count(&status, MPI_DOUBLE, &receiveBufferSize);
  	std::vector<char> receiveBuffer(receiveBufferSize);

  	MPI_Recv(receiveBuffer.data(), receiveBufferSize, MPI_CHAR, source, 0, communicator, MPI_STATUS_IGNORE);

		size_t sizeOfParticleAttributes = sizeof(Attributes);
		char* receiveBufferStart = &receiveBuffer[0];
		char* receiveBufferEnd = &receiveBuffer[receiveBufferSize-1];
  	for (char* i = receiveBufferStart; i != receiveBufferEnd; i += sizeOfParticleAttributes){
			Attributes* attributes = reinterpret_cast<Attributes*>(i);

    	ParticleType particle;
			applyAttributestoParticle(*attributes, particle);
    	receivedParticles.push_back(particle);
  	}
	}
}

MDFlexMPI::MDFlexMPI(int argc, char** argv) : MDFlexSimulation(argc, argv){
  MPI_Init(&argc, &argv);

	int dimensionCount = _configuration->boxMin.value.size();
	double* globalBoxMin = &_configuration->boxMin.value[0];
	double* globalBoxMax = &_configuration->boxMax.value[0];
	_domainDecomposition = std::make_unique<RegularGridDecomposition>(dimensionCount, globalBoxMin, globalBoxMax);

  _particlePropertiesLibrary = std::make_unique<ParticlePropertiesLibraryType>(_configuration->cutoff.value);

	std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();

	for (int i = 0; i < dimensionCount; ++i){
		_configuration->boxMin.value[i] = localBoxMin[i];
		_configuration->boxMax.value[i] = localBoxMax[i];
	}

	_autoPasContainer = std::make_unique<autopas::AutoPas<ParticleType>>(std::cout);	
}

MDFlexMPI::~MDFlexMPI() {
	MPI_Finalize();
}

void MDFlexMPI::run(){

	// @todo: make variable part of MDFlexConfig
	int iterationsPerSuperstep = 10;
	int remainingIterations = _configuration->iterations.value;

	for (int i = 0; i < _configuration->iterations.value; i+=iterationsPerSuperstep){
		executeSuperstep(iterationsPerSuperstep);

	//	remainingIterations -= iterationsPerSuperstep;
	}

	//_simulation->printStatistics(autopasContainer);
}

void MDFlexMPI::executeSuperstep(const int iterationsPerSuperstep){
	updateParticles(iterationsPerSuperstep);

	MPI_Barrier(_domainDecomposition->getCommunicator());

	auto [emigrants, updated] = _autoPasContainer->updateContainer();

	if (updated) {
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

void MDFlexMPI::sendEmigrantsToNeighbours(std::vector<ParticleType> &emigrants){
	std::vector<double> localBoxMin = _domainDecomposition->getLocalBoxMin();
	std::vector<double> localBoxMax = _domainDecomposition->getLocalBoxMax();
	int dimensionCount = localBoxMin.size();

	MPI_Comm communicator = _domainDecomposition->getCommunicator();

	std::vector<ParticleType> immigrants;
	std::vector<ParticleType> migrants;

	std::vector<ParticleType> particlesForPreceedingNeighbour;
	std::vector<ParticleType> particlesForSucceedingNeighbour;

	std::array<MPI_Request, 2> sendRequests;
	std::array<std::vector<char>, 2> sendBuffers;

	for(int i = 0; i < dimensionCount; ++i){
		// @todo: neighbours are dependent on type of domain decomposition.
		// Move neighbour identification to domain decomposition class.
		int preceedingNeighbour = dimensionCount * 2;
		int succeedingNeighbour = dimensionCount * 2 + 1;

		for (auto &particle : emigrants){
			if (particle.getR()[i] - localBoxMin[i] < 0){
				particlesForPreceedingNeighbour.push_back(particle);
			}
			else if (particle.getR()[i] - localBoxMax[i] > 0) {
				particlesForSucceedingNeighbour.push_back(particle);
			}
		}

		sendParticles(particlesForPreceedingNeighbour, preceedingNeighbour, communicator, sendRequests[0],
			sendBuffers[0]);
		sendParticles(particlesForPreceedingNeighbour, succeedingNeighbour, communicator, sendRequests[1],
			sendBuffers[1]);

		receiveParticles(migrants, preceedingNeighbour, communicator);
		receiveParticles(migrants, succeedingNeighbour, communicator);

		MPI_Status sendStatus;
		for (auto request : sendRequests){
			MPI_Wait(&request, &sendStatus);
		}

		particlesForPreceedingNeighbour.clear();
		particlesForSucceedingNeighbour.clear();

		sendBuffers[0].clear();
		sendBuffers[1].clear();

		int nextDimensionIndex = (i + 1) % dimensionCount;
		for (auto &particle : migrants){
			if (particle.getR()[nextDimensionIndex] - localBoxMin[nextDimensionIndex] < 0){
				particlesForPreceedingNeighbour.push_back(particle);
			}
			else if (particle.getR()[nextDimensionIndex] - localBoxMax[nextDimensionIndex] > 0 ){
				particlesForSucceedingNeighbour.push_back(particle);
			}
			else {
				immigrants.push_back(particle);
			}
		}

		migrants.clear();

		sendParticles(particlesForPreceedingNeighbour, preceedingNeighbour, communicator, sendRequests[0], sendBuffers[0]);
		sendParticles(particlesForPreceedingNeighbour, succeedingNeighbour, communicator, sendRequests[1], sendBuffers[1]);

		receiveParticles(immigrants, preceedingNeighbour, communicator);
		receiveParticles(immigrants, succeedingNeighbour, communicator);

		for (auto request : sendRequests){
			MPI_Wait(&request, &sendStatus);
		}

		particlesForPreceedingNeighbour.clear();
		particlesForSucceedingNeighbour.clear();

		sendBuffers[0].clear();
		sendBuffers[1].clear();
	}

	for (auto &particle : immigrants) {
		_autoPasContainer->addParticle(particle);
	}
}
