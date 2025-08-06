#pragma once

#include <math.h>
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"

#include <iostream>

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
            double dist2 = autopas::utils::ArrayMath::dot(dr, dr);
            if (dist2 > cutoffSquared) {
                return;
            }

            double dist = std::sqrt(dist2);
            double adist = alpha * dist;
            double f1 = M_2_SQRTPI * alpha * exp(- adist * adist);
            double f2 = erfc(adist) / dist;

            // does not include coulomb constant
            auto f = dr * i.getQ() * j.getQ() * ((f1 + f2)/dist2);
            //if(i.isOwned())
            //    std::cout << "Particle " << i.getQ() << " short Range F: " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;
            i.addF(f);
            //if(i.isOwned())
            //    std::cout << "Particle " << i.getQ() << " total F: " << i.getF()[0] << ", " << i.getF()[1] << ", " << i.getF()[2] << std::endl;
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