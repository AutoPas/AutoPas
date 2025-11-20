#pragma once
#include "PairwiseFunctor.h"
#include "autopas/utils/kokkos/KokkosSoAType.h"

namespace autopas {

template <class Particle_T, class CRTP_T>
class KokkosFunctor : public PairwiseFunctor<Particle_T, CRTP_T> {
 public:
  using SoAArraysType = typename autopas::utils::kokkos::KokkosSoAType<typename Particle_T::SoAArraysType>;

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
  void KokkosSoAFunctor(SoAArraysType _soa, size_t index) {
    utils::ExceptionHandler::exception("{}::KokkosSoAFunctor: not implemented", this->getName());
  };
};

}  // namespace autopas