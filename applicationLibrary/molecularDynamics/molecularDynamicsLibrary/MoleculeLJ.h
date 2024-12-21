#pragma once

#include "MoleculeLJBase.h"

namespace mdLib {
/**
 * MoleculeLJ with all variables in 32 bit precision
 */
using MoleculeLJFP32 = MoleculeLJBase<float, float, unsigned long>;
/**
 * MoleculeLJ with all variables in 64 bit precision
 */

using MoleculeLJFP64 = MoleculeLJBase<double, double, unsigned long>;
/**
 * MoleculeLJ with mixed precision for calculations and accumulation
 */
using MoleculeLJMP = MoleculeLJBase<float, double, unsigned long>;

#if AUTOPAS_PRECISION_MODE == SPSP
using MoleculeLJ = MoleculeLJFP32;
#elif AUTOPAS_PRECISION_MODE == DPDP
using MoleculeLJ = MoleculeLJFP64;
#else
using MoleculeLJ = MoleculeLJMP;
#endif
}  // namespace mdLib