/**
 * @file Functor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <type_traits>

#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/SoAView.h"
#include "autopas/utils/logging/FLOPLogger.h"

namespace autopas {

/**
 * Newton 3 modes for the Functor.
 */
enum class FunctorN3Modes {
  Newton3Only,
  Newton3Off,
  Both,
};

/**
 * Functor base class.
 * Both an array of structure (AoS) and a structure of array (SoA) are supported
 * to be used with functors.
 * Newton3: A functor does not have to implement both a newton3 and a
 * non-newton3 version. Instead you can specify, which version you use by
 * overriding allowsNonNewton3 resp. allowsNewton3
 *
 * @tparam Particle_T the type of Particle
 * @tparam CRTP_T the actual type of the functor
 */
template <class Particle_T, class CRTP_T>
class Functor {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Make the Implementation type template publicly available.
   */
  using Functor_T = CRTP_T;
  /**
   * Constructor
   * @param cutoff
   */
  explicit Functor(double cutoff) : _cutoff(cutoff){};

  virtual ~Functor() = default;

  /**
   * This function is called at the start of each traversal.
   * Use it for resetting global values or initializing them.
   */
  virtual void initTraversal(){};

  /**
   * This function is called at the end of each traversal.
   * You may accumulate values in this step.
   * @param newton3
   */
  virtual void endTraversal(bool newton3){};

  /**
   * Get attributes needed for computation.
   * @return Attributes needed for computation.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 0> getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 0>{};
  }

  /**
   * Get attributes needed for computation without N3 optimization.
   * @return Attributes needed for computation.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 0> getNeededAttr(std::false_type) {
    return Functor_T::getNeededAttr();
  }

  /**
   * Get attributes computed by this functor.
   * @return Attributes computed by this functor.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 0>{};
  }

  /**
   * Copies the AoS data of the given cell in the given soa.
   *
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added to the SoA with the specified offset.
   * @param skipSoAResize If resizing of the SoA buffers should be skipped or not. If this is called with true, it must
   * be ensured before the call that there is sufficient capacity in the SoA.
   * @note The parameter skipSoAResize is usually set to false, only for VerletListsCellBased Containers it is set to
   * true, since they resize the SoA before the call to SoALoader.
   * @tparam ParticleCell Type of the cell.
   */
  template <class ParticleCell>
  void SoALoader(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset, bool skipSoAResize) {
    SoALoaderImpl(cell, soa, offset, skipSoAResize, std::make_index_sequence<Functor_T::getNeededAttr().size()>{});
  }

  /**
   * Copies the data stored in the soa back into the cell.
   *
   * @param cell Cell where the data should be stored.
   * @param soa  Structure of arrays from where the data is loaded.
   * @param offset Offset within the SoA. The data of the soa should be
   * extracted starting at offset.
   * @tparam ParticleCell Type of the cell.
   */
  template <typename ParticleCell>
  void SoAExtractor(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset) {
    SoAExtractorImpl(cell, soa, offset, std::make_index_sequence<Functor_T::getComputedAttr().size()>{});
  }

  /**
   * Specifies whether the functor is capable of Newton3-like functors.
   * If the functor provides an interface to soa or aos functions that utilize
   * Newton's third law of motion (actio = reactio) to reduce the computational
   * complexity this function should return true. If this is not the case this
   * function should return false.
   * @return true if and only if this functor provides an interface to
   * Newton3-like functions.
   */
  virtual bool allowsNewton3() = 0;

  /**
   * Specifies whether the functor is capable of non-Newton3-like functors.
   * If the functor provides an interface to soa or aos functions that do not
   * utilize Newton's third law of motion (actio = reactio) this function should
   * return true. If this is not the case this function should return false.
   * @return true if and only if this functor provides an interface to functions
   * that do not utilize Newton3.
   */
  virtual bool allowsNonNewton3() = 0;

  /**
   * Specifies whether the functor is capable of using the specified Vectorization Pattern in the SoA functor.
   *
   * @param vecPattern
   * @return whether the functor is capable of using the specified Vectorization Pattern
   */
  virtual bool isVecPatternAllowed(const VectorizationPatternOption::Value vecPattern) = 0;

  /**
   * Specifies whether the functor should be considered for the auto-tuning process.
   * @return true if and only if this functor is relevant for auto-tuning.
   */
  virtual bool isRelevantForTuning() = 0;

  /**
   * Returns name of functor. Intended for use with the iteration logger, to differentiate between calls to
   * computeInteractions using different functors in the logs.
   * @return name of functor.
   */
  virtual std::string getName() = 0;

  /**
   * Getter for the functor's cutoff
   * @return
   */
  [[nodiscard]] double getCutoff() const { return _cutoff; }

  /**
   * Setter for the vectorization pattern to be used
   * @param vecPattern
   */
  virtual void setVecPattern(const VectorizationPatternOption::Value vecPattern) {}

  /**
   * Get the number of FLOPs. Implementation required if FLOPLogger used.
   *
   * If derived class provides no implementation, the FLOPLogger interprets the default numeric_limits<size_t>::max()
   * output as invalid and leaves "Not Implemented" the log.
   * @return number of FLOPs
   */
  [[nodiscard]] virtual size_t getNumFLOPs() const { return std::numeric_limits<size_t>::max(); }

  /**
   * Get the hit rate. Implementation required if FLOPLogger used.
   *
   * If derived class provides no implementation, the FLOPLogger interprets the default NaN output as invalid and
   * leaves "Not Implemented" in the log.
   *
   * @return (number of kernel calls) / (number of distance calculations)
   */
  [[nodiscard]] virtual double getHitRate() const { return std::numeric_limits<double>::quiet_NaN(); }

 private:
  /**
   * Implements loading of SoA buffers.
   * @tparam cell_t Cell type.
   * @tparam I Attribute.
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added
   * @param skipSoAResize If resizing of the SoA buffers should be skipped or not. If this is called with true, it must
   * be ensured before the call that there is sufficient capacity in the SoA.
   * to the SoA with the specified offset.
   */

  template <typename cell_t, std::size_t... I>
  void SoALoaderImpl(cell_t &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset, bool skipSoAResize,
                     std::index_sequence<I...>) {
    if (not skipSoAResize) {
      soa.resizeArrays(offset + cell.size());
    }

    if (cell.isEmpty()) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesn't understand that pointer is used later
    [[maybe_unused]] auto const pointer = std::make_tuple(soa.template begin<Functor_T::getNeededAttr()[I]>()...);

    auto cellIter = cell.begin();
    // load particles in SoAs
    for (size_t i = offset; cellIter != cell.end(); ++cellIter, ++i) {
      /**
       * The following statement writes the values of all attributes defined in neededAttr into the respective position
       * inside the SoA buffer. I represents the index inside neededAttr. The whole expression is folded sizeof...(I)
       * times over the comma operator. E.g. like this (std::index_sequence<I...> = 0, 1):
       * ((std::get<0>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[0]>()),
       * (std::get<1>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[1]>()))
       */
      ((std::get<I>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[I]>()), ...);
    }
  }

  /**
   * Implements extraction of SoA buffers.
   * @tparam cell_t Cell type.
   * @tparam I Attribute.
   * @param cell Cell.
   * @param soa SoA buffer.
   * @param offset Offset
   */
  template <typename cell_t, std::size_t... I>
  void SoAExtractorImpl(cell_t &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset, std::index_sequence<I...>) {
    if (cell.isEmpty()) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesn't understand that pointer is used later
    [[maybe_unused]] auto const pointer = std::make_tuple(soa.template begin<Functor_T::getComputedAttr()[I]>()...);

    auto cellIter = cell.begin();
    // write values in SoAs back to particles
    for (size_t i = offset; cellIter != cell.end(); ++cellIter, ++i) {
      /**
       * The following statement writes the value of all attributes defined in computedAttr back into the particle.
       * I represents the index inside computedAttr.
       * The whole expression is folded sizeof...(I) times over the comma operator. E.g. like this
       * (std::index_sequence<I...> = 0, 1):
       * (cellIter->template set<Functor_T::getComputedAttr()[0]>(std::get<0>(pointer)[i]),
       * cellIter->template set<Functor_T::getComputedAttr()[1]>(std::get<1>(pointer)[i]))
       */
      (cellIter->template set<Functor_T::getComputedAttr()[I]>(std::get<I>(pointer)[i]), ...);
    }
  }

  double _cutoff;
};

}  // namespace autopas