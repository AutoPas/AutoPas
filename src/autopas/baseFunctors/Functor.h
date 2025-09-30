/**
 * @file Functor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <type_traits>

#include "autopas/baseFunctors/FunctorBase.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * Functor base class with CRTP.
 * Both an array of structures (AoS) and a structure of arrays (SoA) are supported to be used with functors.
 * The AoS interface is defined via the pure virtual function AoSFunctor, which must be implemented by the derived
 *
 * @tparam Particle_T the type of Particle
 * @tparam CRTP_T the actual type of the functor
 */
template <class Particle_T, class CRTP_T>
class Functor : public FunctorBase {
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
   * @param name
   */
  explicit Functor(double cutoff, std::string name = "Functor") : FunctorBase(cutoff, name){};

  ~Functor() override = default;

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
   * @tparam ParticleCell_T Type of the cell.
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added to the SoA with the specified offset.
   * @param skipSoAResize If resizing of the SoA buffers should be skipped or not. If this is called with true, it must
   * be ensured before the call that there is sufficient capacity in the SoA.
   * @note The parameter skipSoAResize is usually set to false. Only for VerletListsCellBased Containers it is set to
   * true, since they resize the SoA before the call to SoALoader.
   */
  template <class ParticleCell_T>
  void SoALoader(ParticleCell_T &cell, SoA<SoAArraysType> &soa, size_t offset, bool skipSoAResize) {
    SoALoaderImpl(cell, soa, offset, skipSoAResize, std::make_index_sequence<Functor_T::getNeededAttr().size()>{});
  }

  /**
   * Copies the data stored in the soa back into the cell.
   *
   * @tparam ParticleCell_T Type of the cell.
   * @param cell Cell where the data should be stored.
   * @param soa  Structure of arrays from where the data is loaded.
   * @param offset Offset within the SoA. The data of the soa should be extracted starting at offset.
   */
  template <typename ParticleCell_T>
  void SoAExtractor(ParticleCell_T &cell, SoA<SoAArraysType> &soa, size_t offset) {
    SoAExtractorImpl(cell, soa, offset, std::make_index_sequence<Functor_T::getComputedAttr().size()>{});
  }

 private:
  /**
   * Implements loading of SoA buffers.
   * @tparam ParticleCell_T Cell type.
   * @tparam I Attribute.
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added to the SoA with the specified offset.
   * @param skipSoAResize If resizing of the SoA buffers should be skipped or not. If this is called with true, it must
   * be ensured before the call that there is sufficient capacity in the SoA.
   */
  template <typename ParticleCell_T, std::size_t... I>
  void SoALoaderImpl(ParticleCell_T &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset, bool skipSoAResize,
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
   * @tparam ParticleCell_T Cell type.
   * @tparam I Attribute.
   * @param cell Cell.
   * @param soa SoA buffer.
   * @param offset Offset
   */
  template <typename ParticleCell_T, std::size_t... I>
  void SoAExtractorImpl(ParticleCell_T &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset,
                        std::index_sequence<I...>) {
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
};

}  // namespace autopas