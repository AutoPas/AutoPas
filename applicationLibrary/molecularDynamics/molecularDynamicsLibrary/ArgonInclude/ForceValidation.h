/**
* @file ArgonFunctor.h
* @author I. Angelucci
* @date 02/07/24
*/

#pragma once

#include "RepulsiveTerm.h"
#include "DispersionTerm.h"
#include "DisplacementHandle.h"

namespace autopas::utils::ArrayMath::Argon {

static constexpr std::array<double, 23> A{
    {-0.170806873130E+01,  // 000
     -0.316818997395E+02,  // 001
     -0.571545817498E+05,  // 011
     0.848780677578E+02,   // 111
     0.163923794220E+07,   // 002
     0.380809830366E+02,   // 012
     -0.217403993198E+03,  // 112
     0.244997545538E+03,   // 022
     0.128926029735E+03,   // 122
     0.815601247450E+02,   // 222
     0.409987725022E+02,   // 003
     -0.978512983041E+06,  // 013
     0.104383189893E+07,   // 113
     -0.383820796134E+02,  // 023
     0.143934125087E+03,   // 123
     0.102161665959E+04,   // 033
     -0.569593762549E+02,  // 004
     0.178356269769E+04,   // 014
     0.242202158097E+02,   // 114
     -0.279617357863E+01,  // 024
     -0.324585542907E+02,  // 005
     -0.963264559888E-01,  // 015
     -0.898942588279E+05}  // 006
};

static constexpr std::array<double, 23> alpha{
    {0.428132039316E+00,  // 000
     0.503934786518E+00,  // 001
     0.104706730543E+01,  // 011
     0.456769339560E+00,  // 111
     0.131047310452E+01,  // 002
     0.444052360076E+00,  // 012
     0.480469535570E+00,  // 112
     0.737327026170E+00,  // 022
     0.496177745527E+00,  // 122
     0.424365319847E+00,  // 222
     0.428946186456E+00,  // 003
     0.117979281352E+01,  // 013
     0.119534448663E+01,  // 113
     0.416753172892E+00,  // 023
     0.507114743788E+00,  // 123
     0.764351644551E+00,  // 033
     0.422619330972E+00,  // 004
     0.757543022081E+00,  // 014
     0.482734248672E+00,  // 114
     0.419340374650E+00,  // 024
     0.635761316281E+00,  // 005
     0.375600311119E+00,  // 015
     0.130334333132E+01}  // 006
};

static constexpr std::array<double, 5> Z{
    {0.273486414323E+03,   // 111
     -0.213475877256E+05,  // 112
     0.108226781130E+07,   // 122
     -0.213710093072E+07,  // 222
     0.364515182541E+06}   // 113
};

static constexpr std::array<double, 5> beta{
    {0.211602562917E+02,  // 111
     0.149623190559E+01,  // 112
     0.132161541056E+01,  // 122
     0.208199482789E+01,  // 222
     0.179870559008E+01}  // 113
};

int main() {
  const double ANGSTROM = 10e-10;

  enum ID{I, J, K};

  // equilateral triangle
  std::array<double, 3> positionI{{0, 0, 0}};
  std::array<double, 3> positionJ{{0, 1, 0}};
  std::array<double, 3> positionK{{std::sqrt(3)/2, 0.5, 0}};

  const double R_min = 0.5;
  const double R_max = 6;
  const double step_size = 0.25;
  for (double r = R_min; r <= R_max; r += step_size) {

    positionI *= r;
    positionJ *= r;
    positionK *= r;

    auto displacementHandleIJ = DisplacementHandle(positionI, positionJ, I, J);
    auto displacementHandleJK = DisplacementHandle(positionJ, positionK, J, K);
    auto displacementHandleKI = DisplacementHandle(positionK, positionI, K, I);

    const auto dispersionPotential = U_dispersive(Z, beta, displacementHandleIJ, displacementHandleJK, displacementHandleKI);
    const auto repulsivePotential = U_repulsive(A, alpha, displacementHandleIJ, displacementHandleJK, displacementHandleKI);
    const auto totalPotential = dispersionPotential + repulsivePotential;

    const auto dispersiveForceI = F_dispersive<I>(Z, beta, displacementHandleIJ, displacementHandleJK, displacementHandleKI);
    //const auto repulsiveForceI = F_repulsive<I>(A, alpha, displacementHandleIJ, displacementHandleJK, displacementHandleKI);

  }
}
} // namespace autopas::utils::ArrayMath::Argon
