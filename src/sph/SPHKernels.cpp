//
// Created by seckler on 22.01.18.
//
#include "autopas.h"
#include "SPHKernels.h"

#include <cmath>

const double pi = atan(1.0) * 4.0;
const double kernelSupportRadius = 2.5;

double autopas::sph::W(const std::array<double, 3> dr, const double h) {
    const double H = kernelSupportRadius * h;
    const double s = sqrt(autopas::arrayMath::dot(dr, dr)) / H;  // sqrt(dr * dr) / H;
    const double s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
    const double s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
    double r_value = (s1 * s1 * s1) - 4.0 * (s2 * s2 * s2);  //  pow(s1, 3) - 4.0 * pow(s2, 3);
    //if # of dimension == 3
    r_value *= 16.0 / pi / (H * H * H);
    return r_value;
}

std::array<double, 3> autopas::sph::gradW(const std::array<double, 3> dr, const double h) {
    const double H = kernelSupportRadius * h;
    const double drabs = sqrt(autopas::arrayMath::dot(dr, dr));
    const double s = drabs / H;  // sqrt(dr * dr) / H;
    const double s1 = (1.0 - s < 0) ? 0 : 1.0 - s;
    const double s2 = (0.5 - s < 0) ? 0 : 0.5 - s;
    double r_value = -3.0 * (s1 * s1) + 12.0 * (s2 * s2);  // - 3.0 * pow(s1, 2) + 12.0 * pow(s2, 2);
    //if # of dimension == 3
    r_value *= 16.0 / pi / (H * H * H);
    const double scale = r_value / (drabs * H + 1.0e-6 * h);
    return autopas::arrayMath::mulScalar(dr, scale);  // dr * r_value / (sqrt(dr * dr) * H + 1.0e-6 * h);
}