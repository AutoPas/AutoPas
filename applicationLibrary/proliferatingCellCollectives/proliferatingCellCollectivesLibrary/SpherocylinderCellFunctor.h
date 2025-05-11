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

  // For demonstration, we use length as a proxy for density.
  // In a real implementation, add a density member to SpherocylinderCell.

  inline void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    if (i.isDummy() or j.isDummy()) return;

    auto collisionInfo = i.getCollisionInfo(j);

    if (collisionInfo) {
      auto [overlap, normal, contactPoint] = collisionInfo.value();

      // TODO: friction parameter
      double FRICTION = 15;

      // TODO better stress model
      double stressPenalty = exp(overlap) - exp(-0.5);
      i.setStress(i.getStress() + stressPenalty);
      if (newton3) {
        j.setStress(j.getStress() + stressPenalty);
      }

      // Calculate the force vector based on overlap and normal (use friction parameter)
      std::array<double, 3> forceVector = normal * (FRICTION * std::pow(overlap, 1.5));
      i.setF(i.getF() + forceVector);
      if (newton3) {
        j.setF(j.getF() - forceVector);
      }

      // Calculate the torque based on the contact point and normal
      std::array<double, 3> r1 = contactPoint - i.getR();
      std::array<double, 3> torque1 = cross(r1, forceVector);
      i.setTorque(i.getTorque() + torque1);

      if (newton3) {
        std::array<double, 3> r2 = contactPoint - j.getR();
        std::array<double, 3> torque2 = cross(r2, forceVector * (-1.0));
        j.setTorque(j.getTorque() - torque2);
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