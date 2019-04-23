/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * class of helpers for verlet lists
 * @tparam Particle
 */
template <class Particle>
class VerletListHelpers {
  using floatPrecision = typename Particle::ParticleFloatingPointType;

 public:
  /// AOS verlet list storage
  typedef std::unordered_map<Particle *, std::vector<Particle *>> AoS_verletlist_storage_type;

  /// typedef for soa's of verlet list's linked cells (only id and position needs to be stored)
  typedef typename utils::SoAType<size_t, floatPrecision, floatPrecision, floatPrecision>::Type SoAArraysType;

  /// attributes for soa's of verlet list's linked cells (only id and position needs to be stored)
  enum AttributeNames : int { id, posX, posY, posZ };

  /// typedef for verlet-list particle cell type
  typedef FullParticleCell<Particle, SoAArraysType> VerletListParticleCellType;

  /**
   * This functor can generate verlet lists using the typical pairwise
   * traversal.
   */
  class VerletListGeneratorFunctor : public autopas::Functor<Particle, VerletListParticleCellType, SoAArraysType> {
    typedef VerletListParticleCellType ParticleCell;
    using floatPrecision = typename Particle::ParticleFloatingPointType;

   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoffskin
     */
    VerletListGeneratorFunctor(AoS_verletlist_storage_type &verletListsAoS, floatPrecision cutoffskin)
        : _verletListsAoS(verletListsAoS), _cutoffskinsquared(cutoffskin * cutoffskin) {}

    bool isRelevantForTuning() override { return false; }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      floatPrecision distsquare = ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared)
        // this is thread safe, only if particle i is accessed by only one
        // thread at a time. which is ensured, as particle i resides in a
        // specific cell and each cell is only accessed by one thread at a time
        // (ensured by traversals)
        // also the list is not allowed to be resized!

        _verletListsAoS.at(&i).push_back(&j);
    }

    /**
     * SoAFunctor for verlet list generation. (single cell version)
     * @param soa the soa
     * @param newton3 whether to use newton 3
     */
    void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
      if (soa.getNumParticles() == 0) return;

      auto **const __restrict__ idptr = reinterpret_cast<Particle **const>(soa.template begin<AttributeNames::id>());
      floatPrecision *const __restrict__ xptr = soa.template begin<AttributeNames::posX>();
      floatPrecision *const __restrict__ yptr = soa.template begin<AttributeNames::posY>();
      floatPrecision *const __restrict__ zptr = soa.template begin<AttributeNames::posZ>();

      size_t numPart = soa.getNumParticles();
      for (unsigned int i = 0; i < numPart; ++i) {
        auto &currentList = _verletListsAoS.at(idptr[i]);

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const floatPrecision drx = xptr[i] - xptr[j];
          const floatPrecision dry = yptr[i] - yptr[j];
          const floatPrecision drz = zptr[i] - zptr[j];

          const floatPrecision drx2 = drx * drx;
          const floatPrecision dry2 = dry * dry;
          const floatPrecision drz2 = drz * drz;

          const floatPrecision dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList.push_back(idptr[j]);
            if (not newton3) {
              _verletListsAoS.at(idptr[j]).push_back(idptr[i]);
            }
          }
        }
      }
    }

    /**
     * SoAFunctor for the verlet list generation. (two cell version)
     * @param soa1 soa of first cell
     * @param soa2 soa of second cell
     */
    void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool /*newton3*/) override {
      if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

      auto **const __restrict__ id1ptr = reinterpret_cast<Particle **const>(soa1.template begin<AttributeNames::id>());
      floatPrecision *const __restrict__ x1ptr = soa1.template begin<AttributeNames::posX>();
      floatPrecision *const __restrict__ y1ptr = soa1.template begin<AttributeNames::posY>();
      floatPrecision *const __restrict__ z1ptr = soa1.template begin<AttributeNames::posZ>();

      auto **const __restrict__ id2ptr = reinterpret_cast<Particle **const>(soa2.template begin<AttributeNames::id>());
      floatPrecision *const __restrict__ x2ptr = soa2.template begin<AttributeNames::posX>();
      floatPrecision *const __restrict__ y2ptr = soa2.template begin<AttributeNames::posY>();
      floatPrecision *const __restrict__ z2ptr = soa2.template begin<AttributeNames::posZ>();

      size_t numPart1 = soa1.getNumParticles();
      for (unsigned int i = 0; i < numPart1; ++i) {
        auto &currentList = _verletListsAoS.at(id1ptr[i]);

        size_t numPart2 = soa2.getNumParticles();

        for (unsigned int j = 0; j < numPart2; ++j) {
          const floatPrecision drx = x1ptr[i] - x2ptr[j];
          const floatPrecision dry = y1ptr[i] - y2ptr[j];
          const floatPrecision drz = z1ptr[i] - z2ptr[j];

          const floatPrecision drx2 = drx * drx;
          const floatPrecision dry2 = dry * dry;
          const floatPrecision drz2 = drz * drz;

          const floatPrecision dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList.push_back(id2ptr[j]);
          }
        }
      }
    }

    /**
     * SoALoader for verlet list generator.
     * Only loads IDs and positions
     * @param cell
     * @param soa
     * @param offset
     */
    void SoALoader(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset = 0) override {
      assert(offset == 0);
      soa.resizeArrays(cell.numParticles());

      if (cell.numParticles() == 0) return;

      unsigned long *const __restrict__ idptr = soa.template begin<AttributeNames::id>();
      floatPrecision *const __restrict__ xptr = soa.template begin<AttributeNames::posX>();
      floatPrecision *const __restrict__ yptr = soa.template begin<AttributeNames::posY>();
      floatPrecision *const __restrict__ zptr = soa.template begin<AttributeNames::posZ>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = 0; cellIter.isValid(); ++cellIter, ++i) {
        Particle *pptr = &(*cellIter);
        idptr[i] = reinterpret_cast<std::uintptr_t>(pptr);
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
      }
    }

    /**
     * SoAExtractor for verlet list generation.
     * Currently empty.
     * @param cell
     * @param soa
     * @param offset
     */
    void SoAExtractor(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset = 0) override {
      // nothing yet...
    }

   private:
    AoS_verletlist_storage_type &_verletListsAoS;
    floatPrecision _cutoffskinsquared;
  };

  /**
   * This functor checks the validity of neighborhood lists.
   * If a pair of particles has a distance of less than the cutoff radius it
   * checks whether the pair is represented in the verlet list.
   * If the pair is not present in the list the neigborhood lists are invalid
   * and neighborlistsAreValid()  will return false.
   * @todo: SoA?
   * @tparam ParticleCell
   */
  template <class ParticleCell>
  class VerletListValidityCheckerFunctor : public autopas::Functor<Particle, ParticleCell, SoAArraysType> {
    using floatPrecision = typename Particle::ParticleFloatingPointType;

   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoffsquared
     */
    VerletListValidityCheckerFunctor(AoS_verletlist_storage_type &verletListsAoS, floatPrecision cutoffsquared)
        : _verletListsAoS(verletListsAoS), _cutoffsquared(cutoffsquared), _valid(true) {}

    bool isRelevantForTuning() override { return false; }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      floatPrecision distsquare = ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffsquared) {
        // this is thread safe, we have variables on the stack
        auto found = std::find(_verletListsAoS[&i].begin(), _verletListsAoS[&i].end(), &j);
        if (found == _verletListsAoS[&i].end()) {
          // this is thread safe, as _valid is atomic
          _valid = false;
        }
      }
    }

    /**
     * Returns whether the neighbour list are valid.
     * Call this after performing the pairwise traversal
     * @return
     */
    bool neighborlistsAreValid() { return _valid; }

   private:
    AoS_verletlist_storage_type &_verletListsAoS;
    floatPrecision _cutoffsquared;

    // needs to be thread safe
    std::atomic<bool> _valid;
  };

};  // class VerletListHelpers
}  // namespace autopas
