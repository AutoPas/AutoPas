//
// Created by markus on 9/25/25.
//

#pragma once
#include "autopas/utils/SoA.h"
#include "FunctorWrapperInterface.h"
#include "Functor.h"

template <class Particle_T, class Functor_T>
class FunctorWrapper : public FunctorWrapperInterface<Particle_T> {
public:

  FunctorWrapper(Functor_T f) : functor(std::move(f)) {}

  virtual ~FunctorWrapper() = default;

  void initTraversal() override {functor.initTraversal();}

  void endTraversal(bool newton3) override {functor.endTraversal(newton3);}

  template<class ParticleCell_T>
  void SoALoader(ParticleCell_T &cell, void *soaRaw, size_t offset, bool skipSoAResize) {
    auto &soa = *reinterpret_cast<autopas::SoA<typename Functor_T::SoAArraysType>*>(soaRaw);
    functor.SoALoader(cell, soa, offset, skipSoAResize);
  }

  template<class ParticleCell_T>
  void SoAExtractor(ParticleCell_T &cell, void *soaRaw, size_t offset) {
        auto &soa = *reinterpret_cast<autopas::SoA<typename Functor_T::SoAArraysType>*>(soaRaw);
        functor.SoAExtractor(cell, soa, offset);
   }

  bool allowsNewton3() override {return functor.allowsNewton3();}
  bool allowsNonNewton3() override {return functor.allowsNonNewton3();}

  bool isRelevantForTuning() override {return functor.isRelevantForTuning();}

  std::string getName() override {return functor.getName();}

  [[nodiscard]] double getCutoff() const override {return functor.getCutoff();}

  [[nodiscard]] size_t getNumFlops() const override {return functor.getNumFlops();}

  [[nodiscard]] double getHitRate() const override {return functor.getHitRate();}

  protected:
   Functor_T functor;
};