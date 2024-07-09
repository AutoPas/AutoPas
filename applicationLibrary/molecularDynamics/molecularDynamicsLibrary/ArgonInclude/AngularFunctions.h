/**
 * @file ArgonFunctor.h
 * @author I. Angelucci
 * @date 19/06/24
 */

#pragma once

#include "CosineHandle.h"
#include "DisplacementHandle.h"

namespace autopas::utils::ArrayMath::Argon {

/**
 *
 * @param a first l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @param b second l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @param c third l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @param cosineI angle between the displacements IJ and IK
 * @param cosineJ angle between the displacements JK and JI
 * @param cosineK angle between the displacements KI and KJ
 * @return angular function derived by Bell describing the interaction between dipole/quadrupole/octupole
 */
[[nodiscard]] double AngularTerm(const size_t a, const size_t b, const size_t c, const DisplacementHandle &displacementIJ, const DisplacementHandle &displacementJK,
                                 const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                 const CosineHandle &cosineJ, const CosineHandle &cosineK);

/**
 *
 * @tparam a first l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam b second l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam c third l-value of the angular function, either 1 (dipole), 2 (quadrupole), 3 (octupole)
 * @tparam ID ID of the particle with respect to whose position we are calculating the derivative
 * @param displacementIJ displacement between particle I and particle J
 * @param displacementJK displacement between particle J and particle K
 * @param displacementKI displacement between particle K and particle I
 * @param cosineI angle between the displacements IJ and IK
 * @param cosineJ angle between the displacements JK and JI
 * @param cosineK angle between the displacements KI and KJ
 * @return derivative of angular function derived by Bell describing the interaction between dipole/quadrupole/octupole
 */
template <size_t a, size_t b, size_t c, size_t ID>
[[nodiscard]] double AngularTerm_derive_wrt(const DisplacementHandle &displacementIJ,
                                            const DisplacementHandle &displacementJK,
                                            const DisplacementHandle &displacementKI, const CosineHandle &cosineI,
                                            const CosineHandle &cosineJ, const CosineHandle &cosineK);

}  //  namespace autopas::utils::ArrayMath::Argon