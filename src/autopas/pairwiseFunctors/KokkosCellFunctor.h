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


        template<class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3 = true>
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




        template<class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3>::processCell(
                ParticleCell &cell) {
#ifdef AUTOPAS_KOKKOS
            auto particles = cell._particles;

            if(DataLayout==DataLayoutOption::aos){

                /*
                Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, particles.size()), KOKKOS_LAMBDA(const int i) {
                    for (unsigned int l = i + 1; l < particles.size(); l++) {
                        _functor->AoSFunctorInline(particles[i], particles[l]);
                        _functor->AoSFunctorInline(particles[l], particles[i]);
                    }
                });
    */
                //new implementation - newton 3 off
                if(useNewton3){
                    //no parallelization possible
                    for (unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = i + 1; l < particles.size(); l++) {
                            _functor->AoSFunctorInline(particles[i], particles[l], useNewton3);
                            //_functor->AoSFunctorInline(particles[l], particles[i], useNewton3);
                        }
                    }
                }else{
                    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, particles.size()), KOKKOS_LAMBDA(const unsigned int i) {
                        for (unsigned int l = 0; l < particles.size(); l++) {
                            if(i != l) {
                                _functor->AoSFunctorInline(particles[i], particles[l], useNewton3);
                            }
                        }
                    });
                }//end newton3 off
            }else if(DataLayout == DataLayoutOption::kokkos){
                if(useNewton3){
                    //no parallelization possible
                    for (unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = i + 1; l < particles.size(); l++) {
                            _functor->SoAFunctorInline(cell, cell, i, l,  useNewton3);
                            //_functor->AoSFunctorInline(particles[l], particles[i], useNewton3);
                        }
                    }
                }else{
                    //Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, particles.size())
                    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, particles.size()), KOKKOS_LAMBDA(const unsigned int i) {
                        for (unsigned int l = 0; l < particles.size(); l++) {
                            if(i != l) {
                                _functor->SoAFunctorInline(cell, cell, i, l,  useNewton3);
                            }
                        }
                    });
                }//end newton3 off


            }
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3>::processCellPair(
                ParticleCell &cell1, ParticleCell &cell2) {
#ifdef AUTOPAS_KOKKOS
            auto part0 = cell1._particles;
            auto part1 = cell2._particles;
            /*
            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, part0.size()), KOKKOS_LAMBDA(const int i) {
                for (unsigned int l = 0; l < part1.size(); l++) {
                    _functor->AoSFunctorInline(part0[i], part1[l], useNewton3);
                    _functor->AoSFunctorInline(part1[l], part0[i], useNewton3);
                }
            });*/
            if(DataLayout == DataLayoutOption::aos){


                if(useNewton3){
                    //no parallelization possible
                    for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                        for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                            _functor->AoSFunctorInline(cell1._particles[i], cell2._particles[l], useNewton3);
                        }
                    }
                }else{
                    //Newton 3 off
                    Kokkos::parallel_for(2, KOKKOS_LAMBDA(const int x){
                        if(x == 0){
                            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, part0.size()), KOKKOS_LAMBDA(const unsigned int i) {
                                for (unsigned int l = 0; l < part1.size(); l++) {
                                    _functor->AoSFunctorInline(part0[i], part1[l], useNewton3);

                                }
                            });
                        }else {
                            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, part1.size()), KOKKOS_LAMBDA(const unsigned int i) {
                                for (unsigned int l = 0; l < part0.size(); l++) {
                                    _functor->AoSFunctorInline(part1[i], part0[l], useNewton3);
                                }
                            });
                        }
                    });
                }//end newton3 off
            }else if(DataLayout == DataLayoutOption::kokkos){
                if(useNewton3){
                    //no parallelization possible
                    for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                        for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                            _functor->SoAFunctorInline(cell1, cell2, i, l, useNewton3);
                        }
                    }
                }else{
                    //Newton 3 off
                    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 2), KOKKOS_LAMBDA(const int x){
                        if(x == 0){
                            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, part0.size()), KOKKOS_LAMBDA(const unsigned int i) {
                                for (unsigned int l = 0; l < part1.size(); l++) {
                                    _functor->SoAFunctorInline(cell1, cell2, i, l, useNewton3);
                                }
                            });
                        }else {
                            Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, part1.size()), KOKKOS_LAMBDA(const unsigned int i) {
                                for (unsigned int l = 0; l < part0.size(); l++) {
                                    _functor->SoAFunctorInline(cell2, cell1, i, l, useNewton3);
                                }
                            });
                        }
                    });
                }//end newton3 off

            }
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3>::processCellOneThread(ParticleCell &cell) {
#ifdef AUTOPAS_KOKKOS
            auto particles = cell._particles;
            //std::cout << particles.size() << "\n";
            //always uses newton3
            if(DataLayout == DataLayoutOption::aos){
                if(useNewton3){
                    for (unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = i + 1; l < particles.size(); l++) {
                                _functor->AoSFunctorInline(particles[i], particles[l], true);
                            //_functor->AoSFunctorInline(particles[l], particles[i], useNewton3);
                        }
                    }
                }else{
                    //no parallel for used here
                    for(unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = 0; l < particles.size(); l++) {
                            if(i != l) {
                                _functor->AoSFunctorInline(particles[i], particles[l], useNewton3);
                            }
                        }
                    };
                }
            }else if(DataLayout == DataLayoutOption::kokkos){
                if(useNewton3){
                    for (unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = i + 1; l < particles.size(); l++) {
                            _functor->SoAFunctorInline(cell, cell, i, l, true);
                            //_functor->AoSFunctorInline(particles[l], particles[i], useNewton3);
                        }
                    }
                }else{
                    //no parallel for used here
                    for(unsigned int i = 0; i < particles.size(); i++) {
                        for (unsigned int l = 0; l < particles.size(); l++) {
                            if(i != l) {
                                _functor->SoAFunctorInline(cell, cell, i, l, useNewton3);
                            }
                        }
                    };
                }
            }
#endif
        }

        template<class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3>
        void KokkosCellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3>::processCellPairOneThread(ParticleCell &cell1,
                                                                                                  ParticleCell &cell2) {
#ifdef AUTOPAS_KOKKOS
            //Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i) {
            if(DataLayout == DataLayoutOption::aos) {
                for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                    for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                        if (useNewton3) {
                            _functor->AoSFunctorInline(cell1._particles[i], cell2._particles[l], true);
                        } else {
                            _functor->AoSFunctorInline(cell1._particles[i], cell2._particles[l], useNewton3);
                            _functor->AoSFunctorInline(cell2._particles[l], cell1._particles[i], useNewton3);
                        }
                        //_functor->AoSFunctorInline(cell2._particles[l], cell1._particles[i]);
                    }
                }
            }else if(DataLayout == DataLayoutOption::kokkos){

                if (useNewton3) {
                    for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                        for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                                _functor->SoAFunctorInline(cell1, cell2, i, l, true);
                        }
                    }
                }else{
                    for (unsigned int i = 0; i < cell1._particles.size(); i++) {
                        for (unsigned int l = 0; l < cell2._particles.size(); l++) {
                            _functor->SoAFunctorInline(cell1, cell2, i, l, useNewton3);
                            _functor->SoAFunctorInline(cell2, cell1, l, i, useNewton3);
                        }
                    }
                }//end newton3
            } else {

            }
            //});
#endif


        }
    }
}