/**
 * @file NearFieldFunctor.h
 * @date 30.09.19
 * @author Joachim Marin
 */

#pragma once

#include "Math3D.h"
#include "FmmParticle.h"
#include "autopas/pairwiseFunctors/Functor.h"

class NearFieldFunctor : public autopas::Functor<FmmParticle, autopas::FullParticleCell<FmmParticle>> {

 public:
  explicit NearFieldFunctor(typename FmmParticle::ParticleFloatingPointType cutoff) : Functor(cutoff) {

  }

  void AoSFunctor(FmmParticle &i, FmmParticle &j, bool newton3) override {
    auto distVec = Math3D::subtract(i.getR(), j.getR());
    auto spherical = Math3D::toSpherical(distVec);
    auto dist = spherical[0];
    i.resultFMM += j.charge / dist;
  }

  bool allowsNewton3() override {return false;}
  bool allowsNonNewton3() override {return true;}
  bool isRelevantForTuning() override {return false;}

};

