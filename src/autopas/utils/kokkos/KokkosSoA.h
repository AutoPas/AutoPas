#pragma once

#include "autopas/utils/SoAStorage.h"

namespace autopas::utils::kokkos {

template <class SoAArraysType>
class KokkosSoA {
public:
  KokkosSoA() = default;
  KokkosSoA(const KokkosSoA& soa) = default;


private:

  SoAStorage<SoAArraysType> _soaStorage;
};


}