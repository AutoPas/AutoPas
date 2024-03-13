/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "OpenMPConfigurator.h"

namespace autopas {
/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <long>
 */
    unsigned long openMPDefaultChunkSize = 1;

/**
 * OpenMP configurator default constructor.
 */
    autopas::OpenMPConfigurator::OpenMPConfigurator() = default;

/**
 * OpenMP configurator constructor.
 * @param s chunk size used in OpenMP's loop scheduling
 */
    autopas::OpenMPConfigurator::OpenMPConfigurator(unsigned long s) {
        _chunkSize = s;
    }

/**
 * OpenMP chunk size getter.
 * @return the current OpenMP chunk size
 */
    unsigned long autopas::OpenMPConfigurator::getChunkSize() const {
        return _chunkSize;
    }

/**
 * OpenMP chunk size setter.
 */
    void autopas::OpenMPConfigurator::setChunkSize(unsigned long s) {
        _chunkSize = s;
    }
}