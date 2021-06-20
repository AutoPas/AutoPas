/**
 * @file Checkpoint.h
 * @author N. Fottner
 * @date 10/9/19
 */
#pragma once

#include <fstream>
#include <sstream>
#include <string>

/**
 * In this namespace the initialization of an AutoPas container from vtk checkpoint files is implemented.
 */
namespace Checkpoint {

// anonymous namespace for namespace-private helper functions
namespace {
/**
 * Reads the next numberOfParticles lines from file. The lines are expected to contain xyz data.
 * @tparam dataType type of the data to be read. (Can be scalars or containers that support [])
 * @tparam size number of entries per dataType.
 * @param file
 * @param numberOfParticles
 * @return Vector of read data.
 */
template <class dataType, int size>
std::vector<dataType> readPayload(std::ifstream &file, size_t numberOfParticles) {
  std::vector<dataType> data(numberOfParticles);
  // loop over every line (=particle)
  for (size_t i = 0; i < numberOfParticles; ++i) {
    // loop over line (=coordinates)
    if constexpr (size == 1) {
      file >> data[i];
    } else {
      for (size_t j = 0; j < size; ++j) {
        file >> data[i][j];
      }
    }
  }
  return data;
}
}  // namespace

/**
 * Searches the file word by word and sets the file accessor directly behind the first found position.
 * @param file
 * @param word
 */
void findWord(std::ifstream &file, const std::string &word) {
  std::string currentWord;
  while (currentWord != word) {
    file >> currentWord;
  }
}

/**
 * Reads the Data of all particles from a Vtk file and adds them into the AutoPas container.
 * @param autopas
 * @param vtkFilename
 */
template <class AutoPasTemplate>
void loadParticles(AutoPasTemplate &autopas, const std::string &vtkFilename) {
  std::ifstream infile(vtkFilename);
  size_t numParticles;
  std::string dataType;

  findWord(infile, "POINTS");
  infile >> numParticles >> dataType;

  auto positions = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  // next payload block is always preceded by the datatype
  findWord(infile, dataType);
  auto velocities = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  findWord(infile, dataType);
  auto forces = readPayload<std::array<double, 3>, 3>(infile, numParticles);
  findWord(infile, "default");
  auto typeID = readPayload<size_t, 1>(infile, numParticles);

  // creating Particles from checkpoint:
  for (auto i = 0ul; i < numParticles; ++i) {
    typename AutoPasTemplate::Particle_t p;
    using namespace autopas::utils::ArrayUtils;
    p.setR(static_cast_array<typename AutoPasTemplate::Particle_t::ParticleSoAFloatPrecision>(positions[i]));
    p.setV(static_cast_array<typename AutoPasTemplate::Particle_t::ParticleSoAFloatPrecision>(velocities[i]));
    p.setF(static_cast_array<typename AutoPasTemplate::Particle_t::ParticleSoAFloatPrecision>(forces[i]));
    p.setTypeId(typeID[i]);
    p.setID(i);
    autopas.addParticle(p);
  }
}
}  // namespace Checkpoint
