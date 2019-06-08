/**
 * @file DataLayoutOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

namespace autopas {
/**
 * Possible choices for the particle data layout.
 */
enum DataLayoutOption { aos, soa, cuda };

/**
 * Provides a way to iterate over the possible choices of data layouts.
 */
static const std::set<DataLayoutOption> allDataLayoutOptions = {
    DataLayoutOption::aos,
    DataLayoutOption::soa,
#if defined(AUTOPAS_CUDA)
    DataLayoutOption::cuda,
#endif
};

}  // namespace autopas
