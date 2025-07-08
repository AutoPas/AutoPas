#pragma once

#include <math.h>
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas{
    template <class Particle_T>
    class P3M_shortRangeFunctor : public PairwiseFunctor<Particle_T, P3M_shortRangeFunctor<Particle_T>>{

        private:
        double alpha;
        double cutoffSquared;

        P3M_shortRangeFunctor(): PairwiseFunctor<Particle_T, P3M_shortRangeFunctor<Particle_T>>(0.0), alpha(0.0), cutoffSquared(0.0){}

        public:
        P3M_shortRangeFunctor(double alpha, double cutoff): PairwiseFunctor<Particle_T, P3M_shortRangeFunctor<Particle_T>>(cutoff), alpha(alpha), cutoffSquared(cutoff*cutoff){}

        void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
            using namespace autopas::utils::ArrayMath::literals;
            
            if (i.isDummy() or j.isDummy()) {
                return;
            }
            auto dr = i.getR() - j.getR();
            double dist = autopas::utils::ArrayMath::dot(dr, dr);
            if (dist > cutoffSquared) {
                return;
            }

            double adist = alpha * dist;
            double f1 = M_2_SQRTPI * alpha * exp(- adist * adist);
            double f2 = erfc(adist) / dist;

            auto f = dr * i.getQ() * j.getQ() * ((f1 + f2)/dist);
            i.addF(f);
            if (newton3) {
                // only if we use newton 3 here, we want to
                j.subF(f);
            }
        }

        bool allowsNewton3(){
            return true;
        }
        
        bool allowsNonNewton3(){
            return true;
        }

        bool isRelevantForTuning(){
            return false;
        }

         std::string getName(){
            return "shortRangeP3M";
         }
    };
}