/**
 * @file KokkosFlopCounterFunctor.h
 * @author M. Geitner
 * @date 22.07.19
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

    template<class Particle, class ParticleCell>
    class KokkosFlopCounterFunctor : public Functor<Particle, ParticleCell> {
        typedef typename Particle::SoAArraysType SoAArraysType;

    public:
        bool isRelevantForTuning() override { return false; }

        /**
         * constructor of FlopCounterFunctor
         * @param cutoffRadius the cutoff radius
         */
        explicit KokkosFlopCounterFunctor<Particle, ParticleCell>(double cutoffRadius)
                : autopas::Functor<Particle, ParticleCell>(),
                  _cutoffSquare(cutoffRadius * cutoffRadius),
                  _distanceCalculations(0ul),
                  _kernelCalls(0ul) {}


                  //INLINE
        void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
          //@TODO change to kokkos components
          auto dr = ArrayMath::sub(i.getR(), j.getR());
          double dr2 = ArrayMath::dot(dr, dr);


          ++_distanceCalculations;

          if (dr2 <= _cutoffSquare) ++_kernelCalls;



        }


        /**
       * Empty SoAExtractor.
       * Nothing to be done yet.
       */
        AUTOPAS_FUNCTOR_SOAEXTRACTOR(, , ,)

        /**
         * get the hit rate of the pair-wise interaction, i.e. the ratio of the number
         * of kernel calls compared to the number of distance calculations
         * @return the hit rate
         */
        double getHitRate() { return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations); }

        /**
         * get the total number of flops
         * @param numFlopsPerKernelCall
         * @return
         */
        double getFlops(unsigned long numFlopsPerKernelCall) const {
          const double distFlops = numFlopsPerDistanceCalculation * static_cast<double>(_distanceCalculations);
          const double kernFlops = numFlopsPerKernelCall * static_cast<double>(_kernelCalls);
          return distFlops + kernFlops;
        }

        /**
         * get the number of calculated distance operations
         * @return
         */
        unsigned long getDistanceCalculations() const { return _distanceCalculations; }

        /**
         * get the number of kernel calls, i.e. the number of pairs of particles with
         * a distance not larger than the cutoff
         * @return
         */
        unsigned long getKernelCalls() const { return _kernelCalls; }

        /**
         * number of flops for one distance calculation.
         * 3 sub + 3 square + 2 add
         */
        static constexpr double numFlopsPerDistanceCalculation = 8.0;


        void initTraversal() override {

        }

        void endTraversal(bool newton3) override {

        }

        void SoAFunctor(SoA <SoAArraysType> &soa, bool newton3) override {

        }

        void SoAFunctor(SoA <SoAArraysType> &soa,
                        const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList,
                        size_t iFrom, size_t iTo, bool newton3) override {

        }

        void SoAFunctor(SoA <SoAArraysType> &soa1, SoA <SoAArraysType> &soa2, bool newton3) override {

        }

        void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) override {

        }

        void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                         CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) override {

        }

        void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                             CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {

        }

        void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                                CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) override {

        }

        void SoALoader(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset) override {

        }

        bool allowsNewton3() override {
          return false;
        }

        bool allowsNonNewton3() override {
          return true;
        }

    private:
        double _cutoffSquare;
        unsigned long _distanceCalculations, _kernelCalls;

    };


}