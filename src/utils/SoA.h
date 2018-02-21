#ifndef AUTOPAS_SOA_H
#define AUTOPAS_SOA_H

#include <map>
#include <vector>
#include "AlignedAllocator.h"

namespace autopas {

class SoA {
 public:
  /**
   * @brief Default constructor.
   */
  SoA() = default;

  /**
   * @brief Copy constructor.
   * @param soa SoA to copy.
   */
  SoA(const SoA &soa) = default;

  /**
   * @brief Destructor.
   */
  ~SoA() {
    for (auto p : arrays) {
      delete p.second;
    }
  }

  /**
   * @brief Creates an aligned vector for every given attribute.
   * @param attributes Vector of Attributes that shall be stored.
   */
  template<std::size_t arraySize>
  void initArrays(std::array<int, arraySize> attributes) {
    for (auto a : attributes) {
      // TODO: assert attribute value does not already exit
      arrays.insert(
          make_pair(a, new std::vector<double, AlignedAllocator<double>>));
    }
  }

  /**
   * @brief Pushes a given value to the desired attribute array.
   * @param attribute Index of array to push to.
   * @param value Value to push.
   */
  // TODO: for maximum convenience make vectors pushable
  void push(const int attribute, const double value) {
    arrays[attribute]->push_back(value);
  }

  /**
   * @brief Reads from all given attribute arrays at position `particleId`.
   * @tparam ArrayLength length of the returned array. Should be equal
   * attributes.size().
   * @param attributes Attributes to read from.
   * @param particleId Position to read from.
   * @return Array of attributes ordered by given attribute order.
   */
  template<std::size_t arrayLength>
  std::array<double, arrayLength> read(std::array<int, arrayLength> attributes,
                                       unsigned int particleId) {
    std::array<double, arrayLength> retArray;
    int i = 0;
    for (auto a : attributes) {
      retArray[i++] = arrays[a]->at(particleId);
    }
    return retArray;
  }

  /**
   * @brief Reads the value of a given attribute of a given particle.
   * @param attributes Attribute to read from.
   * @param particleId Position to read from.
   * @return Attribute value.
   */
  double read(int attribute, unsigned int particleId) {
    return arrays[attribute]->at(particleId);
  }

  double *begin(int attribute) {
    return &(arrays[attribute]->front());
  }

  /**
   * @brief Writes / updates values of attributes for a specific particle.
   * @tparam ArrayLength length of the attributes and value array.
   * @param attributes Array of attributes to update.
   * @param particleId Particle to update.
   * @param value New value.
   */
  template<int arrayLength>
  void write(std::array<int, arrayLength> attributes, unsigned int particleId,
             std::array<double, arrayLength> value) {
    int i = 0;
    for (auto a : attributes) {
      arrays[a]->at(particleId) = value[i++];
    }
  }

  template<int arrayLength>
  void add(std::array<int, arrayLength> attributes, unsigned int particleId,
           std::array<double, arrayLength> value) {
    int i = 0;
    for (auto a : attributes) {
      arrays[a]->at(particleId) += value[i++];
    }
  }

  template<int arrayLength>
  void sub(std::array<int, arrayLength> attributes, unsigned int particleId,
           std::array<double, arrayLength> value) {
    int i = 0;
    for (auto a : attributes) {
      arrays[a]->at(particleId) -= value[i++];
    }
  }

  /**
   * @brief Returns the number of particles.
   *
   * This function only checks the size of the first array since it is assumed
   * that the user manages the arrays responsibly.
   *
   * @return Number of particles.
   */
  size_t getNumParticles() const {
    if (arrays.empty()) {
      return 0;
    }
    return arrays.begin()->second->size();
  }

 private:
  // TODO: variable precision (two maps?, user defined via initArrays?)
  // TODO: maybe fix number of attributes via template?
  std::map<int, std::vector<double, AlignedAllocator<double>> *> arrays;
};
}  // namespace autopas
#endif  // AUTOPAS_SOA_H
