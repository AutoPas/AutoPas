//
// Created by maximilian on 12.08.19.
//

#pragma once

namespace autopas{
    namespace utils{
    template <class Functor, DataLayoutOption dataLayout>
    class KokkosDataLayoutConverter{
    public:
        /**
         * Constructor
         * @tparam Functor Functor Type
         * @param functor responsible for the conversion
         */
        KokkosDataLayoutConverter(Functor *functor) : _functor(functor) {}

        /**
         * loads the target dataLayout in a cell
         * @tparam ParticleCell Cell type
         * @param cell to load the data in
         */
        template <class ParticleCell>
        void loadDataLayout(ParticleCell &cell) {
            _functor->kokkosExtractor(cell);
        }

        /**
         * converts the dataLayout to aos
         * @tparam ParticleCell Cell type
         * @param cell to load the data in
         */
        template <class ParticleCell>
        void storeDataLayout(ParticleCell &cell) {
            _functor->kokkosLoader(cell);
        }

    private:
        /**
         *  Functor to convert cells
         */
        Functor *_functor;
    };
    }//namespace utils
}//namespace autopas