/**
 * @file VerletListHelpers.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>
namespace autopas {

/**
 * class of helpers for verlet lists
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class VerletListHelpers {
 public:
  /// AOS verlet list storage
  typedef std::unordered_map<Particle *, std::vector<Particle *>> AoS_verletlist_storage_type;

  /**
   * This functor can generate verlet lists using the typical pairwise
   * traversal.
   * @todo: SoA?
   */
  class VerletListGeneratorFunctor : public autopas::Functor<Particle, ParticleCell> {
   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoffskinsquared
     */
    VerletListGeneratorFunctor(AoS_verletlist_storage_type &verletListsAoS, double cutoffskinsquared)
        : _verletListsAoS(verletListsAoS), _cutoffskinsquared(cutoffskinsquared) {}

    void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      double distsquare = ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared)
        // this is thread safe, only if particle i is accessed by only one
        // thread at a time. which is ensured, as particle i resides in a
        // specific cell and each cell is only accessed by one thread at a time
        // (ensured by traversals)
        // also the list is not allowed to be resized!

        _verletListsAoS.at(&i).push_back(&j);
    }

    virtual void SoAFunctor(SoA<Particle> &soa, bool /*newton3*/ = true) override {
      if (soa.getNumParticles() == 0) return;

      Particle **const __restrict__ idptr =
          reinterpret_cast<Particle **const>(soa.template begin<Particle::AttributeNames::id>());
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

      size_t numPart = soa.getNumParticles();
      for (unsigned int i = 0; i < numPart; ++i) {
        auto &currentList = _verletListsAoS.at(idptr[i]);

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const double drx = xptr[i] - xptr[j];
          const double dry = yptr[i] - yptr[j];
          const double drz = zptr[i] - zptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList.push_back(idptr[j]);
          }
        }
      }
    }

    virtual void SoAFunctor(SoA<Particle> &soa1, SoA<Particle> &soa2, bool /*newton3*/ = true) override {
      if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

      Particle **const __restrict__ id1ptr =
          reinterpret_cast<Particle **const>(soa1.template begin<Particle::AttributeNames::id>());
      double *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();

      Particle **const __restrict__ id2ptr =
          reinterpret_cast<Particle **const>(soa2.template begin<Particle::AttributeNames::id>());
      double *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

      size_t numPart1 = soa1.getNumParticles();
      for (unsigned int i = 0; i < numPart1; ++i) {
        auto &currentList = _verletListsAoS.at(id1ptr[i]);

        size_t numPart2 = soa2.getNumParticles();

        for (unsigned int j = 0; j < numPart2; ++j) {
          const double drx = x1ptr[i] - x2ptr[j];
          const double dry = y1ptr[i] - y2ptr[j];
          const double drz = z1ptr[i] - z2ptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList.push_back(id2ptr[j]);
          }
        }
      }
    }

    virtual void SoALoader(ParticleCell &cell, SoA<Particle> &soa, size_t offset = 0) override {
      assert(offset == 0);
      soa.resizeArrays(cell.numParticles());

      if (cell.numParticles() == 0) return;

      unsigned long *const __restrict__ idptr = soa.template begin<Particle::AttributeNames::id>();
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = 0; cellIter.isValid(); ++cellIter, ++i) {
        Particle *pptr = &(*cellIter);
        idptr[i] = reinterpret_cast<unsigned long &>(pptr);
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
      }
    }

    virtual void SoAExtractor(ParticleCell &cell, SoA<Particle> &soa, size_t offset = 0) override {
      // nothing yet...
    }

   private:
    AoS_verletlist_storage_type &_verletListsAoS;
    double _cutoffskinsquared;
  };

  /**
   * This functor checks the validity of neighborhood lists.
   * If a pair of particles has a distance of less than the cutoff radius it
   * checks whether the pair is represented in the verlet list.
   * If the pair is not present in the list the neigborhood lists are invalid
   * and neighborlistsAreValid()  will return false.
   * @todo: SoA?
   */
  class VerletListValidityCheckerFunctor : public autopas::Functor<Particle, ParticleCell> {
   public:
    /**
     * Constructor
     * @param verletListsAoS
     * @param cutoffsquared
     */
    VerletListValidityCheckerFunctor(AoS_verletlist_storage_type &verletListsAoS, double cutoffsquared)
        : _verletListsAoS(verletListsAoS), _cutoffsquared(cutoffsquared), _valid(true) {}

    void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      double distsquare = ArrayMath::dot(dist, dist);
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
    double _cutoffsquared;

    // needs to be thread safe
    std::atomic<bool> _valid;
  };

};  // class VerletListHelpers
}  // namespace autopas