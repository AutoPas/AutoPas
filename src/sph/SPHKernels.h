//
// Created by seckler on 22.01.18.
//

#ifndef AUTOPAS_SPHKERNELS_H
#define AUTOPAS_SPHKERNELS_H

#include <array>

namespace autopas {
    namespace sph {
        double W(const std::array<double, 3> dr, const double h);
        unsigned long getFlopsW();
        std::array<double, 3> gradW(const std::array<double, 3> dr, const double h);
    }  // namespace autopas
}  // namespace autopas
#endif //AUTOPAS_SPHKERNELS_H
