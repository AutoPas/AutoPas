/**
 * @file RegularGrid.cpp
 * @author J. Körner
 * @date 19.04.2021
 */
#include "RegularGrid.h"

#include <algorithm>
#include <functional>
#include <list>
#include <math.h>
#include <numeric>

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

	void calculatePrimeFactors(unsigned int number, std::list<unsigned int>& oPrimeFactors){
		while (number%2 == 0)
		{
			oPrimeFactors.push_back(2);
			number = number / 2;
		}

		for (unsigned int i = 3; i <= number; i = i+2)
		{
			while (number%i == 0)
			{
				oPrimeFactors.push_back(i);
				number = number / i;
			}
		}
	}

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
		attributes.oldForceX = particle.getOldF()[0];
		attributes.oldForceY = particle.getOldF()[1];
		attributes.oldForceZ = particle.getOldF()[2];
	}

	void applyAttributestoParticle(Attributes &attributes, ParticleType &particle){
		particle.setR({attributes.positionX, attributes.positionY, attributes.positionZ});
		particle.setV({attributes.velocityX, attributes.velocityY, attributes.velocityZ});
		particle.setF({attributes.forceX, attributes.forceY, attributes.forceZ});
		particle.setOldF({attributes.oldForceX, attributes.oldForceY, attributes.oldForceZ});
		particle.setID(attributes.id);
		particle.setTypeId(attributes.typeId);
	}

	const char* serializeAttributes(ParticleType &particle, Attributes &attributes){
 		return reinterpret_cast<const char*>(&(particle.getR()[0])); 
	}
}

RegularGrid::RegularGrid(int argc, char** argv, const int &dimensionCount,
	const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax){
	_dimensionCount = dimensionCount;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &_subdomainCount);

	initializeDecomposition();

	initializeMPICommunicator();

	initializeLocalDomain();

	initializeGlobalBox(globalBoxMin, globalBoxMax);

	initializeLocalBox();

	initializeNeighbourIds();
}

RegularGrid::~RegularGrid(){
	MPI_Finalize();
}

void RegularGrid::update(){
	updateLocalBox();
}

void RegularGrid::initializeDecomposition(){
	std::list<unsigned int> primeFactors;
	calculatePrimeFactors(_subdomainCount, primeFactors);

	while (primeFactors.size() > _dimensionCount)
	{
		primeFactors.sort();
		auto firstElement = primeFactors.front();
		primeFactors.pop_front();
		primeFactors.front() *= firstElement; 
	}

	_decomposition.resize(_dimensionCount);

	for (auto& dimensionSize : _decomposition)
	{
		if (primeFactors.size() > 0) {
			dimensionSize = primeFactors.front();
			primeFactors.pop_front();
		}
		else {
			dimensionSize = 1;
		}
	}
}

void RegularGrid::initializeMPICommunicator(){
	std::vector<int> periods(_dimensionCount, 1);
  MPI_Cart_create(MPI_COMM_WORLD, _dimensionCount, _decomposition.data(), periods.data(), true, &_communicator);
 	MPI_Comm_rank(_communicator, &_domainIndex);
}

void RegularGrid::initializeLocalDomain(){
	_domainId.resize(_dimensionCount);
  MPI_Comm_rank(_communicator, &_domainIndex);

	std::vector<int> periods(_dimensionCount, 1);
	MPI_Cart_get(_communicator, _dimensionCount, _decomposition.data(), periods.data(), _domainId.data());
}

void RegularGrid::initializeLocalBox(){
	_localBoxMin.resize(_dimensionCount);
	_localBoxMax.resize(_dimensionCount);
	updateLocalBox();
}

void RegularGrid::initializeNeighbourIds(){
	_neighbourDomainIndices.resize(_dimensionCount * 2);

	for (int i = 0; i < _dimensionCount; ++i){
		auto neighbourIndex = i * 2;
		auto preceedingNeighbourId = _domainId;
		preceedingNeighbourId[i] = (--preceedingNeighbourId[i] + _decomposition[i])%_decomposition[i];
		_neighbourDomainIndices[neighbourIndex] = convertIdToIndex(preceedingNeighbourId);

		++neighbourIndex;
		auto succeedingNeighbourId = _domainId;
		succeedingNeighbourId[i] = (++succeedingNeighbourId[i] + _decomposition[i])%_decomposition[i];
		_neighbourDomainIndices[neighbourIndex] = convertIdToIndex(succeedingNeighbourId);
	}
}

int RegularGrid::convertIdToIndex(const std::vector<int> &domainId){
	int neighbourDomainIndex = 0;

	for (int i = 0; i < _dimensionCount; ++i){
		int accumulatedTail = 1;

		if (i < _decomposition.size()-1) {
			accumulatedTail = std::accumulate(_decomposition.begin()+i+1, _decomposition.end(),
				1, std::multiplies<int>());
		}

		neighbourDomainIndex += accumulatedTail * domainId[i];
	}

	return neighbourDomainIndex;
}

void RegularGrid::updateLocalBox(){
  for (int i = 0; i < _dimensionCount; ++i) {
		double localBoxWidth = (_globalBoxMax[i] - _globalBoxMin[i]) / static_cast<double>(_decomposition[i]);

    _localBoxMin[i] = _domainId[i] * localBoxWidth + _globalBoxMin[i];
    _localBoxMax[i] = (_domainId[i] + 1) * localBoxWidth + _globalBoxMin[i];

    if (_domainId[i] == 0) {
      _localBoxMin[i] = _globalBoxMin[i];
    } else if (_domainId[i] == _decomposition[i] - 1) {
      _localBoxMax[i] = _globalBoxMax[i];
    }
  }
}

void RegularGrid::initializeGlobalBox(const std::vector<double> &globalBoxMin,
	const std::vector<double> &globalBoxMax){
	_globalBoxMin.resize(_dimensionCount);
	_globalBoxMax.resize(_dimensionCount);
	for (int i = 0; i < _dimensionCount; ++i){
		_globalBoxMin[i] = globalBoxMin[i];
		_globalBoxMax[i] = globalBoxMax[i];
	}
}

void RegularGrid::exchangeHaloParticles(SharedAutoPasContainer &autoPasContainer) {
	// @todo: create configuration parameter for halo width
	double haloWidth;

	int dimensionCount = _localBoxMin.size();
	int neighbourCount = dimensionCount * 2;

	std::vector<ParticleType> particlesForLeftNeighbour;
	std::vector<ParticleType> particlesForRightNeighbour;
	std::vector<ParticleType> haloParticles;

	double particlePosition;
	
	// @todo: these are named in a misinforming manner. They do not refer to the current domain's halos,
	// but to the halo boundaries of the respective neighbour
	double leftNeighbourHaloEnd, rightNeighbourHaloStart;
	double bottomNeighbourHaloEnd, topNeighbourHaloStart;

	int leftNeighbour, rightNeighbour;

	int nextDimensionIndex;

	for(int i = 0; i < dimensionCount; ++i){
		leftNeighbour = (i * 2) % neighbourCount;
		rightNeighbour = (i * 2 + 1) % neighbourCount;

		haloWidth = (_localBoxMax[i] - _localBoxMin[i]) / 20.0;
		leftNeighbourHaloEnd = _localBoxMin[i] + haloWidth;
		rightNeighbourHaloStart = _localBoxMax[i] - haloWidth;

		nextDimensionIndex = (i + 1) % dimensionCount;
		// @todo: remove as soon es the configuration parameter for the halo width has been defined
		haloWidth = (_localBoxMax[nextDimensionIndex] - _localBoxMin[nextDimensionIndex]) / 20.0;
		for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
			particlePosition = particle->getR()[i];
			if (particlePosition < leftNeighbourHaloEnd){
				particlesForLeftNeighbour.push_back(*particle);
			}
			else if (particlePosition > rightNeighbourHaloStart) {
				particlesForRightNeighbour.push_back(*particle);
			}
		}
		sendParticles(particlesForLeftNeighbour, leftNeighbour);
		sendParticles(particlesForLeftNeighbour, rightNeighbour);

		receiveParticles(haloParticles, leftNeighbour);
		receiveParticles(haloParticles, rightNeighbour);

		waitForSendRequests();

		particlesForLeftNeighbour.clear();
		particlesForRightNeighbour.clear();

		nextDimensionIndex = (i + 1) % dimensionCount;

		leftNeighbourHaloEnd = _localBoxMin[nextDimensionIndex] + haloWidth;
		rightNeighbourHaloStart = _localBoxMax[nextDimensionIndex] - haloWidth;

		leftNeighbour = (leftNeighbour + 2) % neighbourCount;
		rightNeighbour = (rightNeighbour + 3) % neighbourCount;

		for (auto &particle : haloParticles){
			particlePosition = particle.getR()[nextDimensionIndex];
			if (particlePosition < leftNeighbourHaloEnd){
				particlesForLeftNeighbour.push_back(particle);
			}
			else if (particlePosition > rightNeighbourHaloStart){
				particlesForRightNeighbour.push_back(particle);
			}
		}

		sendParticles(particlesForLeftNeighbour, leftNeighbour);
		sendParticles(particlesForLeftNeighbour, rightNeighbour);

		receiveParticles(haloParticles, leftNeighbour);
		receiveParticles(haloParticles, rightNeighbour);

		waitForSendRequests();

		particlesForLeftNeighbour.clear();
		particlesForRightNeighbour.clear();
	}

	for (auto &particle : haloParticles) {
		autoPasContainer->addOrUpdateHaloParticle(particle);
	}
}

void RegularGrid::exchangeMigratingParticles(SharedAutoPasContainer &autoPasContainer) {
	auto [emigrants, updated] = autoPasContainer->updateContainer();
	int dimensionCount = _localBoxMin.size();

	std::vector<ParticleType> immigrants;
	std::vector<ParticleType> migrants;

	std::vector<ParticleType> particlesForLeftNeighbour;
	std::vector<ParticleType> particlesForRightNeighbour;

	int nextDimensionIndex;
	int leftNeighbour, rightNeighbour;
	int neighbourCount = dimensionCount * 2;
	double particlePosition;

	for(int i = 0; i < dimensionCount; ++i){
		leftNeighbour = (i * 2) % neighbourCount;
		rightNeighbour = (i * 2 + 1) % neighbourCount;

		for (auto &particle : emigrants){
			particlePosition = particle.getR()[i];
			if (particlePosition - _localBoxMin[i] < 0){
				particlesForLeftNeighbour.push_back(particle);
			}
			else if (particlePosition - _localBoxMax[i] > 0) {
				particlesForRightNeighbour.push_back(particle);
			}
		}

		sendParticles(particlesForLeftNeighbour, leftNeighbour);
		sendParticles(particlesForLeftNeighbour, rightNeighbour);

		receiveParticles(migrants, leftNeighbour);
		receiveParticles(migrants, rightNeighbour);

		waitForSendRequests();

		particlesForLeftNeighbour.clear();
		particlesForRightNeighbour.clear();

		leftNeighbour = (leftNeighbour + 2) % neighbourCount;
		rightNeighbour = (rightNeighbour + 3) % neighbourCount;

		nextDimensionIndex = (i + 1) % dimensionCount;

		for (auto &particle : migrants){
			particlePosition = particle.getR()[nextDimensionIndex];
			if (particlePosition - _localBoxMin[nextDimensionIndex] < 0){
				particlesForLeftNeighbour.push_back(particle);
			}
			else if (particlePosition - _localBoxMax[nextDimensionIndex] > 0 ){
				particlesForRightNeighbour.push_back(particle);
			}
			else {
				immigrants.push_back(particle);
			}
		}

		migrants.clear();

		sendParticles(particlesForLeftNeighbour, leftNeighbour);
		sendParticles(particlesForLeftNeighbour, rightNeighbour);

		receiveParticles(immigrants, leftNeighbour);
		receiveParticles(immigrants, rightNeighbour);

		waitForSendRequests();

		particlesForLeftNeighbour.clear();
		particlesForRightNeighbour.clear();
	}

	for (auto &particle : immigrants) {
		autoPasContainer->addParticle(particle);
	}
}

void RegularGrid::sendParticles(std::vector<ParticleType> &particles, int &receiver) {
	std::vector<char> buffer;
 	for (auto &particle : particles) {
		Attributes attributes;
		transferParticleAttributes(particle, attributes);

		const char* attributesStart = serializeAttributes(particle, attributes);
		const char* attributesEnd = attributesStart + sizeof(Attributes) - 1;

   	buffer.insert(std::end(buffer), attributesStart, attributesEnd);
 	}
	sendDataToNeighbour(buffer, receiver);
}

void RegularGrid::receiveParticles(std::vector<ParticleType> &receivedParticles, int &source) {
 	std::vector<char> receiveBuffer;

	receiveDataFromNeighbour(source, receiveBuffer);

	size_t sizeOfParticleAttributes = sizeof(Attributes);
	char* start = &receiveBuffer[0];
	char* end = &receiveBuffer[receiveBuffer.size() - 1];
 	for (char* i = start; i != end; i += sizeOfParticleAttributes){
		Attributes* attributes = reinterpret_cast<Attributes*>(i);

   	ParticleType particle;
		applyAttributestoParticle(*attributes, particle);
   	receivedParticles.push_back(particle);
 	}
}

void RegularGrid::sendDataToNeighbour(std::vector<char> sendBuffer, const int &neighbour){
		_sendBuffers.push_back(sendBuffer);

		MPI_Request sendRequest;
		_sendRequests.push_back(sendRequest);

  	MPI_Isend(sendBuffer.data(), sendBuffer.size(), MPI_CHAR, neighbour, 0, _communicator, &sendRequest);
}

void RegularGrid::receiveDataFromNeighbour(const int &neighbour, std::vector<char> &receiveBuffer){
  	MPI_Status status;
  	MPI_Probe(neighbour, 0, _communicator, &status);

  	int receiveBufferSize;
  	MPI_Get_count(&status, MPI_DOUBLE, &receiveBufferSize);
		receiveBuffer.resize(receiveBufferSize);

  	MPI_Recv(receiveBuffer.data(), receiveBufferSize, MPI_CHAR, neighbour, 0, _communicator, MPI_STATUS_IGNORE);
}

void RegularGrid::waitForSendRequests(){
		MPI_Status sendStatus;
		for (auto request : _sendRequests){
			MPI_Wait(&request, &sendStatus);
		}
		_sendRequests.clear();
		_sendBuffers.clear();
}

void RegularGrid::synchronizeDomains(){
	MPI_Barrier(_communicator);
}
