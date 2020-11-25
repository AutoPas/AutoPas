/**
 * @file namespaces.h
 * In this file all namespaces of autopas should be listed and documented.
 *
 * @author seckler
 * @date 07.06.18
 */

#pragma once

/**
 * This is the main namespace of AutoPas.
 * Everything defined within AutoPas should be in this namespace.
 */
namespace autopas {

/**
 * This namespace is used for implementation specifics.
 * If you are a developer of AutoPas you might want to take a look inside here.
 */
namespace internal {}

/**
 * This namespace is used for memory profiling functions.
 */
namespace memoryProfiler {}

/**
 * Algorithms for creating a smooth function through a series of points.
 */
namespace smoothing {}

/**
 * In this namespace all functionality of the Smoothed Particle Hydrodynamics module of autopas is present.
 * This mainly includes kernels and particles.
 */
namespace sph {}

/**
 * In this namespace some helper classes and functions can be found used inside of AutoPas.
 * These classes reside mostly in the utils directory. However, most commonly used function inside the utils
 * directory might not be in the utils namespace.
 */
namespace utils {

/**
 * Namespace to handle mathematical operations of std::array.
 */
namespace ArrayMath {}

/**
 * In this namespace some helper functions for std::array can be found.
 */
namespace ArrayUtils {}

/**
 * Some functions to parse enums from (input-) strings.
 */
namespace StringUtils {}

/**
 * Namespace to handle the conversion between one dimensional and three dimensional indices.
 * The running index is x.
 */
namespace ThreeDimensionalMapping {}

/**
 * In this namespace some helper functions for std::tuple can be found.
 */
namespace TupleUtils {}

}  // namespace utils

/**
 * Namespace that contains the explicitly defined options of AutoPas.
 */
namespace options {}

}  // namespace autopas
