//
// Created by seckler on 22.01.18.
//

#ifndef AUTOPAS_SPHKERNELS_H
#define AUTOPAS_SPHKERNELS_H

#include <array>

double W(const std::array<double, 3> dr, const double h) ;

std::array<double, 3> gradW(const std::array<double, 3> dr, const double h);


#endif //AUTOPAS_SPHKERNELS_H
