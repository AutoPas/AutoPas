#pragma once
#include "PairwiseFunctor.h"
#include "autopas/utils/kokkos/KokkosSoA.h"
#include "autopas/utils/kokkos/KokkosSoAType.h"

namespace autopas {

template <class Particle_T, class CRTP_T>
class KokkosFunctor : public PairwiseFunctor<Particle_T, CRTP_T> {
 public:
  using SoAArraysType = Particle_T::SoAArraysType;

  explicit KokkosFunctor(double cutoff) : autopas::PairwiseFunctor<Particle_T, CRTP_T>(cutoff) {}

  void SoAFunctorSingle(SoAView<typename PairwiseFunctor<Particle_T, CRTP_T>::SoAArraysType> soa,
                        bool newton3) override {}

  void SoAFunctorVerlet(SoAView<typename PairwiseFunctor<Particle_T, CRTP_T>::SoAArraysType> soa,
                        const size_t indexFirst, const std::vector<size_t, AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override {}
  void SoAFunctorPair(SoAView<typename PairwiseFunctor<Particle_T, CRTP_T>::SoAArraysType> soa1,
                      SoAView<typename PairwiseFunctor<Particle_T, CRTP_T>::SoAArraysType> soa2,
                      bool newton3) override {}

  KOKKOS_INLINE_FUNCTION
  void KokkosSoAFunctor(auto &team, const autopas::utils::kokkos::KokkosSoA<SoAArraysType> &_soa, auto &block,
                        size_t b1_start, size_t b1_end) const {
    utils::ExceptionHandler::exception("{}::KokkosSoAFunctor: not implemented", this->getName());
  };

  KOKKOS_FUNCTION void KokkosSoAFunctorPairwise(auto &team, const utils::kokkos::KokkosSoA<SoAArraysType> &_soa,
                                                auto &block1, size_t i, size_t b1_start, size_t b2_start,
                                                size_t b2_end) const {
    utils::ExceptionHandler::exception("{}::KokkosSoAFunctor: not implemented", this->getName());
  }

  KOKKOS_FUNCTION void load_particle(auto& soa, Particle_T &particle, size_t index) const {
    utils::ExceptionHandler::exception("{}::KokkosSoAFunctor: not implemented", this->getName());
  }

  KOKKOS_FUNCTION void store_particle(auto &soa, size_t index) const {
    utils::ExceptionHandler::exception("{}::KokkosSoAFunctor: not implemented", this->getName());
  }
};

}  // namespace autopas