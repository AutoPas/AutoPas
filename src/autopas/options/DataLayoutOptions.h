/**
 * @file DataLayoutOptions.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

namespace autopas {
/**
 * Possible choices for the particle data layout.
 */
enum DataLayoutOption { aos, soa, cuda };

/**
 * Provides a way to iterate over the possible choices of data layouts.
 */
static const std::vector<DataLayoutOption> allDataLayoutOptions = {
    DataLayoutOption::aos,
    DataLayoutOption::soa,
    DataLayoutOption::cuda,
};

}  // namespace autopas
