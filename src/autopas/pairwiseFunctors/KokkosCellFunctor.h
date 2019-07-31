/**
 *
 * @author M. Geitner
 * @file KokkosCellFunctor.h
 * @date 15.07.19
 *
 */
#pragma once

namespace autopas {
    namespace internal {


        template<class Particle, class ParticleCell, class ParticleFunctor>
        class KokkosCellFunctor {
        public:
            /**
             * The constructor of CellFunctor
             * @param f the particlefunctor which should be used for the interaction.
             */
            explicit KokkosCellFunctor(ParticleFunctor *f) : _functor(f) {}


            /**
             * kokkos functor
             */
            ParticleFunctor *_functor;


            /**
           * process the interactions inside one cell
           * @param cell all pairwise interactions of particles inside this cell are
           * calculated
           */
            void processCell(ParticleCell &cell);

            void processCellOneThread(ParticleCell &cell);

            /**
             * process the interactions between the particles of cell1 with particles of
             * cell2.
             * @param cell1
             * @param cell2
             */
            void processCellPair(ParticleCell &cell1, ParticleCell &cell2);

            void processCellPairOneThread(ParticleCell &cell1, ParticleCell &cell2);
        };


        template<class Particle, class ParticleCell, class ParticleFunctor>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor>::processCell(
                ParticleCell &cell) {
#ifdef AUTOPAS_KOKKOS
            auto particles = cell._particles;
            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, particles.size())/*particles.size()*/, KOKKOS_LAMBDA(const int i) {
                for (unsigned int l = i + 1; l < particles.size(); l++) {
                    _functor->AoSFunctorInline(particles[i], particles[l]);
                    _functor->AoSFunctorInline(particles[l], particles[i]);
                }
            });
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor>::processCellPair(
                ParticleCell &cell1, ParticleCell &cell2) {
#ifdef AUTOPAS_KOKKOS
            auto part0 = cell1._particles;
            auto part1 = cell2._particles;

            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, part0.size()), KOKKOS_LAMBDA(const int i) {
                for (unsigned int l = 0; l < part1.size(); l++) {
                    _functor->AoSFunctorInline(part0[i], part1[l]);
                    _functor->AoSFunctorInline(part1[l], part0[i]);
                }
            });
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor>::processCellOneThread(ParticleCell &cell) {
#ifdef AUTOPAS_KOKKOS
            auto particles = cell._particles;
            //std::cout << particles.size() << "\n";
            Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i) {
                for (unsigned int i = 0; i < particles.size(); i++) {
                    for (unsigned int l = i + 1; l < particles.size(); l++) {
                        _functor->AoSFunctorInline(particles[i], particles[l]);
                        _functor->AoSFunctorInline(particles[l], particles[i]);
                    }
                }
            });
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor>::processCellPairOneThread(ParticleCell &cell1,
                                                                                                  ParticleCell &cell2) {
#ifdef AUTOPAS_KOKKOS
            //Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i) {
                for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                    for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                        _functor->AoSFunctorInline(cell1._particles[i], cell2._particles[l]);
                        _functor->AoSFunctorInline(cell2._particles[l], cell1._particles[i]);
                    }
                }
            //});
#endif


        }
    }
}