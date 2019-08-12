/**
 * @file KokkosTypes.h
 * @author M. Geitner
 * @date 24.06.19
 *
 */

#pragma once

#ifdef AUTOPAS_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace autopas {

/**
 * Change value to change precision
 */
#ifndef KOKKOS_PRECISION
#define KOKKOS_PRECISION 2
#endif
#if KOKKOS_PRECISION == 1
typedef float KOKKOS_FLOAT;
#else
typedef double KOKKOS_FLOAT;
#endif

/**
 * Change value to modify dimensions, currently only support KOKKOS_DIM = 3
 */
#ifndef KOKKOS_DIM
#define KOKKOS_DIM 3
#endif

#ifdef AUTOPAS_KOKKOS
    typedef Kokkos::View<KOKKOS_FLOAT *, Kokkos::DefaultExecutionSpace> FloatVectorType;
    typedef Kokkos::View<KOKKOS_FLOAT *, Kokkos::DefaultExecutionSpace> FloatVectorType;
    typedef Kokkos::View<KOKKOS_FLOAT ***, Kokkos::DefaultExecutionSpace> FloatMatrix3Type;
    typedef Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>  range_policy;
    //typedef Kokkos::View<KOKKOS_FLOAT *> FloatVectorType;
#endif
}  // namespace autopas