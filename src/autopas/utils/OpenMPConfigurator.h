/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

namespace autopas {
/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <long>
 */
extern unsigned long openMPDefaultChunkSize;

/**
 * This class provides configurable parameters for OpenMP.
 */
class OpenMPConfigurator {
    private:
        unsigned long _chunkSize = openMPDefaultChunkSize;

    public:
        /**
         * OpenMP configurator default constructor.
         */
        OpenMPConfigurator();

        /**
         * OpenMP configurator constructor.
         * @param s chunk size used in OpenMP's loop scheduling
         */
        [[maybe_unused]] explicit OpenMPConfigurator(unsigned long s);

        /**
         * OpenMP chunk size getter.
         * @return the current OpenMP chunk size
         */
        [[nodiscard]] unsigned long getChunkSize() const;

        /**
         * OpenMP chunk size setter.
         */
        [[maybe_unused]] void setChunkSize(unsigned long s);
};
}
