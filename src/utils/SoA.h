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
  void initArrays(std::vector<int> attributes) {
    for (auto a : attributes) {
      arrays.insert(
          make_pair(a, new std::vector<double, AlignedAllocator<double>>));
    }
  }

  /**
   * @brief Pushes a given value to the desired attribute array.
   * @param attribute Index of array to push to.
   * @param value Value to push.
   */
  void push(const int attribute, const double value) {
    arrays.at(attribute)->push_back(value);
  }

  /**
   * @brief Reads from all given attribute arrays at position i.
   * @tparam ArrayLength length of the returned array. Should be equal
   * attributes.size().
   * @param attributes Attributes to read from.
   * @param particleId Position to read from.
   * @return Array of attributes ordered by given attribute order.
   */
  template<int arrayLength>
  std::array<double, arrayLength> readParticle(std::vector<int> attributes,
                                               unsigned int particleId) {
    std::array<double, arrayLength> retArray;
    int i = 0;
    for (auto a : attributes) {
      retArray[i++] = arrays.at(a)->at(particleId);
    }
    return retArray;
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
