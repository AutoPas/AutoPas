/**
 * @file SpherocylinderCellFunctor.h
 * @date 11/05/2025
 * @author Manuel Lerchner
 */

#pragma once

#include <cmath>
#include <string>

#include "SpherocylinderCell.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"

namespace pccLib {

using autopas::utils::ArrayMath::operator+=;
using autopas::utils::ArrayMath::operator-=;
using autopas::utils::ArrayMath::operator*;
using autopas::utils::ArrayMath::cross;
using autopas::utils::ArrayMath::operator-;
using autopas::utils::ArrayMath::operator+;

class SpherocylinderCellFunctor : public autopas::PairwiseFunctor<SpherocylinderCell, SpherocylinderCellFunctor> {
 public:
  using Particle = SpherocylinderCell;
  using SoAArraysType = typename Particle::SoAArraysType;
  using AttributeNames = Particle::AttributeNames;

  SpherocylinderCellFunctor() : autopas::PairwiseFunctor<Particle, SpherocylinderCellFunctor>(0.) {}

  std::string getName() override { return "SpherocylinderCellFunctor"; }
  bool isRelevantForTuning() override { return true; }
  bool allowsNewton3() override { return true; }
  bool allowsNonNewton3() override { return true; }

  inline void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    if (i.isDummy() or j.isDummy()) return;

    auto collisionInfo = i.getCollisionInfo(j);

    if (collisionInfo) {
      auto [overlap, normal, qci, qcj] = collisionInfo.value();

      double stressPenalty = 5 * (exp(overlap) - 1);

      double E = 4e6;
      double rc = 1;
      double tc = 0.5;
      double drag = 200;

      std::array<double, 3> forceVector = normal * ((E * tc / (drag)) * std::pow(overlap, 1.5));

      double Li = (i.getLength() / rc);
      std::array<double, 3> torquei = cross(qci - i.getR(), forceVector);
      i.addF(forceVector * (1 / Li));
      i.addTorque(torquei * (12 / std::pow(Li, 3)));
      i.setStress(i.getStress() + stressPenalty);

      if (newton3) {
        double Lj = (j.getLength() / rc);
        std::array<double, 3> torquej = cross(qcj - j.getR(), forceVector);
        j.subF(forceVector * (1 / Lj));
        j.addTorque(torquej * (12 / std::pow(Lj, 3)));
        j.setStress(j.getStress() + stressPenalty);
      }
    }
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) override {
    throw std::runtime_error("SoAFunctorSingle not implemented");
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      bool newton3) override {
    throw std::runtime_error("SoAFunctorPair not implemented");
  }

  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override {
    throw std::runtime_error("SoAFunctorVerlet not implemented");
  }
};

}  // namespace pccLib