/**
* @file KryptonExtendedATMFunctor.h
* @author muehlhaeusser
* @date 08.08.2024
*/

#pragma once

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "LUT3B.h"

namespace mdLib {

/**
* A functor to handle the interactions between three krypton atoms using the extended AxilrodTeller potential.
* This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
* scheme.
* The functor follows the potential described by Jäger et al. in: https://doi.org/10.1063/1.4943959
*
* @note  All calculations assume units to be in Angström and Kelvin.
*
* @tparam Particle The type of particle.
* @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
* @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
* @tparam countFLOPs counts FLOPs and hitrate
*/
    template<class Particle, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool useLUT = false, bool useLUTGlobal = false,
            bool calculateGlobals = false, bool countFLOPs = false>
    class KryptonExtendedATMFunctor
            : public autopas::TriwiseFunctor<Particle,
                    KryptonExtendedATMFunctor<Particle, useNewton3, useLUT, useLUTGlobal, calculateGlobals, countFLOPs>> {
        /**
         * Structure of the SoAs defined by the particle.
         */
        using SoAArraysType = typename Particle::SoAArraysType;

        /**
         * Precision of SoA entries.
         */
        using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

    public:
        /**
         * Deleted default constructor
         */
        KryptonExtendedATMFunctor() = delete;

        /**
         * Constructor
         * @param cutoff
         */
        explicit KryptonExtendedATMFunctor(double cutoff, LUT3B *lut = nullptr)
                : autopas::TriwiseFunctor<Particle,
                KryptonExtendedATMFunctor<Particle, useNewton3, useLUT, useLUTGlobal, calculateGlobals, countFLOPs>>(
                cutoff),
                  _cutoffSquared{cutoff * cutoff},
                  _potentialEnergySum{0.},
                  _virialSum{0., 0., 0.},
                  _aosThreadDataGlobals(),
                  _postProcessed{false} {
            if constexpr (calculateGlobals) {
                _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
            }
            if constexpr (countFLOPs) {
                _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
            }
            if (useLUT) {
                _lut = lut;
            }

        }


        std::string getName() final { return "KryptonExtendedATMFunctor"; }

        bool isRelevantForTuning() final { return true; }

        bool allowsNewton3() final {
            return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
        }

        bool allowsNonNewton3() final {
            return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
        }


        /**
          * Added for consistency, used to fill LUT. Sets the particle properties constants for this functor.
          *
          * This is only necessary if no particlePropertiesLibrary is used.
          *
          *
         */
        void setParticleProperties() {
            //TODO add a way to set NU

            //added for lut by me
            if (useLUT) {

//                _lut->fill<decltype(*this)>(*this, _cutoffSquared, useLUTGlobal);

            }

        }

        void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) override {
            using namespace autopas::utils::ArrayMath::literals;

          newton3 = allowsNewton3();

            if (i.isDummy() or j.isDummy() or k.isDummy()) {
                return;
            }

            const auto threadnum = autopas::autopas_get_thread_num();

            if constexpr (countFLOPs) {
                ++_aosThreadDataFLOPs[threadnum].numDistCalls;
            }

            const auto displacementIJ = j.getR() - i.getR();
            const auto displacementJK = k.getR() - j.getR();
            const auto displacementKI = i.getR() - k.getR();

             double distSquaredIJ = autopas::utils::ArrayMath::dot(displacementIJ, displacementIJ);
             double distSquaredJK = autopas::utils::ArrayMath::dot(displacementJK, displacementJK);
             double distSquaredKI = autopas::utils::ArrayMath::dot(displacementKI, displacementKI);

            double allDistsTripled;


            // Check cutoff for every distance
            if (distSquaredIJ > _cutoffSquared or distSquaredJK > _cutoffSquared or distSquaredKI > _cutoffSquared) {
                return;
            }

            //declare & init so useLUT can use it too
            std::array<double, 3> forceI = {0.,0.,0.};
            std::array<double, 3> forceJ = {0.,0.,0.};
            std::array<double, 3> forceK = {0.,0.,0.};
            double cosines=0;
            double expTerm= 0;
            double sum = 0.0;

            double devIJ , devJK , devKI = 0;
            if(useLUT){


             // std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> res = _lut->retrieveValues(*this, distSquaredIJ, distSquaredKI, distSquaredJK);
              std::pair<const std::array<double, 3>, std::array<u_int8_t, 3>> res = _lut->retrieveValues(*this, distSquaredIJ, distSquaredJK, distSquaredKI);
//              std::pair<const std::array<double, 4>, std::array<u_int8_t, 3>> res = _lut->retrieveValues(*this, distSquaredIJ, distSquaredKI, distSquaredJK);
                auto factors = res.first;
                auto order = res.second;

                //Calculate pre-factor
//                 double distIJ = std::sqrt(distSquaredIJ);
//                 double distJK = std::sqrt(distSquaredJK);
//                 double distKI = std::sqrt(distSquaredKI);
//                 double KIcosIJ = (distSquaredIJ + distSquaredKI - distSquaredJK) / (2 * distIJ * distKI);
//                 double IJcosJK = (distSquaredIJ + distSquaredJK - distSquaredKI) / (2 * distIJ * distJK);
//                 double JKcosKI = (distSquaredJK + distSquaredKI - distSquaredIJ) / (2 * distJK * distKI);
//
//                 double sign_factor = 1 + 3*(KIcosIJ * IJcosJK * JKcosKI);
//
//
//                 // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
//                 const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
//                 const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
//                 const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;
//                 const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
//                 const double numerator = numKI * numJK * numIJ;
//                 double f =  1+  (3. / 8.) * numerator / allDistsSquared;
//
// //                auto A = factors[3];
//                  devIJ = factors[order[0]];
//                  devKI = factors[order[1]];
//                  devJK = factors[order[2]];



//                auto factorForceJDirectionIJ = sign_factor * factors[order[0]];
//                if(factorForceJDirectionIJ<0 ){
//                    factorForceJDirectionIJ = factorForceJDirectionIJ *-1;
//                }
//                auto factorForceJDirectionKI =  sign_factor*  factors[order[1]];;
//                if(factorForceJDirectionKI>0 ){
//                    factorForceJDirectionKI = factorForceJDirectionKI *-1;
//                }
//                auto factorForceJDirectionJK = sign_factor *  factors[order[2]];;
//                if(factorForceJDirectionJK<0 ){
//                    factorForceJDirectionJK = factorForceJDirectionJK *-1;
//                }

//                distSquaredIJ = 0.5;
//                distSquaredKI= 1;
//                distSquaredJK= 0.4;
//
//
//                //testing if order matters
                auto res1 = _lut->getLUTValuesKrypton(distSquaredJK, distSquaredKI, distSquaredIJ);
                auto res2 = _lut->getLUTValuesKrypton(distSquaredJK, distSquaredIJ, distSquaredKI);
                auto res3 = _lut->getLUTValuesKrypton(distSquaredKI, distSquaredIJ, distSquaredJK);
                auto res4 = _lut->getLUTValuesKrypton(distSquaredKI, distSquaredJK, distSquaredIJ);
                auto res5 = _lut->getLUTValuesKrypton(distSquaredIJ, distSquaredJK, distSquaredKI);
                auto res6 = _lut->getLUTValuesKrypton(distSquaredIJ, distSquaredKI, distSquaredJK);


                // Assembling the forces
//                const auto forceIDirectionIJ = displacementIJ * (factorForceJDirectionIJ);
//                const auto forceIDirecationKI = displacementKI * (factorForceJDirectionKI);
//
//                forceI = (forceIDirectionIJ + forceIDirecationKI) * (-1.0);
//
//                //TODO put back
//                i.addF(forceI);
//
//                forceJ = forceI;
//                forceK = forceI;



//new LUT 4 force

                // const auto displacementIJ = j.getR() - i.getR();
                // const auto displacementJK = k.getR() - j.getR();
                // const auto displacementKI = i.getR() - k.getR();
                //
                // const double cosinesGradientIJ =
                //         (3. / 4.) *
                //         ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
                //


// auto forceDirectionIJ = -1. * sign_factor * devIJ * (displacementIJ  );
// auto forceDirectionIJ_2 = (-1. * f) * devIJ * (displacementIJ / distIJ);
// auto forceDirectionKI = -1. * f * devKI * (displacementKI );
// auto forceDirectionKI_2 = (-1. * f) * devKI * (displacementKI / distKI );
// auto forceDirectionJK = -1. * sign_factor * devKI * (displacementJK );
// auto forceDirectionJK_2 = (-1. * f )* devJK * (displacementJK /distJK );

//Assemble force I
              const auto forceIDirectionIJ = displacementIJ* (factors[order[0]]);
              const auto forceIDirectionKI = displacementKI * (factors[order[2]]);
forceI = forceIDirectionIJ - forceIDirectionKI;


// auto forceI_2 = -1. * sign_factor * (devIJ * distSquaredIJ - devKI *distSquaredKI);
//forceI= -1.* (forceDirectionIJ_2 - forceDirectionKI_2);




//auto dcosI_di = (1 /distIJ * distKI) * (distSquaredKI-distSquaredIJ)-
// auto dRi_dri = displacementIJ / distIJ;
// auto dRk_drk = displacementKI / distKI;
// auto force_4 = -1. * (sign_factor *((devIJ* dRi_dri) + (devKI * dRk_drk)));

i.addF(forceI);
                if(newton3){
                    // auto dRi_dri = displacementIJ / distIJ;
                    // auto dRj_drj = displacementJK / distJK;
                    // forceJ = -1. * ( -1. *forceDirectionIJ_2 + forceDirectionJK_2);
                  auto const forceJDirectionIJ = displacementIJ * (factors[order[0]]);
                  auto const forceJDirectionJK = displacementJK * (factors[order[1]]);
                  forceJ=  (forceJDirectionJK - forceJDirectionIJ);

                    j.addF(forceJ);

//                    forceK = -1. * (sign_factor *((devKI* dRk_drk) + (devJK * dRj_drj)));
                    // forceK = -1. * ((-1. * forceDirectionJK_2) + forceDirectionKI_2);
                  auto const forceKDirectionKI = displacementKI * (factors[order[2]]);
                  auto const forceKDirectionJK = displacementJK * (factors[order[1]]);
                  forceK=  (forceI + forceJ) * (-1.0);
                    k.addF(forceK);

                    // Assembling the forces
//                    factorForceJDirectionIJ = -1* factorForceJDirectionIJ;
//
//                    const auto forceJDirectionIJ = displacementIJ * (factorForceJDirectionIJ);
//                    const auto forceJDirectionJK = displacementJK * (factorForceJDirectionJK);
//                    forceJ = (forceJDirectionIJ + forceJDirectionJK) * (-1.0);
//
//                    //TODO
//                    j.addF(forceJ);
//
//                    // Using newton's third law for the force on particle k
//                    forceK = (forceI + forceJ) * (-1.0);
//
//                    //TODO
//                    k.addF(forceK);

                }
                std::cout <<"In LUT should not be here" ;
            }//end of uselut
            else {

                // Actual distances
                const double distIJ = std::sqrt(distSquaredIJ);
                const double distJK = std::sqrt(distSquaredJK);
                const double distKI = std::sqrt(distSquaredKI);

                // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
                const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
                const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
                const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

                const double numerator = numKI * numJK * numIJ;

                const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
                const double allDists = distIJ * distJK * distKI;
                  allDistsTripled = allDistsSquared * allDists;

                // Gradient factors of 1. / (rrr)^3
                const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
                const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

                // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
                cosines = (3. / 8.) * numerator / allDistsSquared;
                const double cosinesGradientIJ =
                        (3. / 4.) *
                        ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
                const double cosinesGradientKI =
                        (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) /
                                     allDistsSquared);

                // Gradient factors corresponding to the normal ATM term
                const auto fullATMGradientIJ =
                        _nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
                const auto fullATMGradientKI =
                        _nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

                expTerm = std::exp(-_alpha * (distIJ + distJK + distKI));

                // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
                std::array<double, 6> sumFactors{};

                for (auto n = 0; n < sumFactors.size(); n++) {
                    sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
                    sum += sumFactors[n];
                }

                // Gradient factor of the sum in ij-direction
                double ijSum = 0.0;
                for (auto n = 0; n < sumFactors.size(); n++) {
                    ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha);
                }

                // Gradient factor of the sum in ki-direction
                double kiSum = 0.0;
                for (auto n = 0; n < sumFactors.size(); n++) {
                    kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha);
                }

                // Total gradient factors for the exponential term times the cosines term
                const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
                const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);

                auto testIJ = expTerm * (devIJ + cosinesGradientIJ * sum);
                auto testKI = expTerm * (devKI + cosinesGradientKI * sum);

                // Assembling the forces
                auto factorIDdirectionIJ = fullATMGradientIJ + fullExpGradientIJ;
                auto factorIDdirectionKI =fullATMGradientKI + fullExpGradientKI;

                const auto forceIDirectionIJ = displacementIJ * (fullATMGradientIJ + fullExpGradientIJ);
                const auto forceIDirecationKI = displacementKI * (fullATMGradientKI + fullExpGradientKI);

                 forceI = (forceIDirectionIJ + forceIDirecationKI) * (-1.0);

                i.addF(forceI);

                 forceJ = forceI;
                 forceK = forceI;

            if (newton3) {
                // Calculate all components for jk-direction
                const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
                const double cosinesGradientJK =
                        (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
                const auto fullATMGradientJK =
                        _nu * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);

                double jkSum = 0.0;
                for (auto n = 0; n < sumFactors.size(); n++) {
                    jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
                }
                const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);



                //delete later
                auto forceNEGJDirectionIJ = (-fullATMGradientIJ - fullExpGradientIJ);
                auto forceNEGJDirectionJK=(fullATMGradientJK + fullExpGradientJK);



                // Assembling the forces
                const auto forceJDirectionIJ = displacementIJ * (-fullATMGradientIJ - fullExpGradientIJ);
                const auto forceJDirectionJK = displacementJK * (fullATMGradientJK + fullExpGradientJK);

                forceJ = (forceJDirectionIJ + forceJDirectionJK) * (-1.0);
                j.addF(forceJ);

                // Using newton's third law for the force on particle k
                forceK = (forceI + forceJ) * (-1.0);
                k.addF(forceK);
            }
            }//End of no_lut

            if constexpr (countFLOPs) {
                if (newton3) {
                    ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
                } else {
                    ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
                }
            }

            if constexpr (calculateGlobals) {
                // Add 3 * potential energy to every owned particle of the interaction.
                // Division to the correct value is handled in endTraversal().
                const double potentialEnergy = (1.0 + cosines) * (_nu / allDistsTripled + expTerm * sum);

                // Virial is calculated as f_i * r_i
                // see Thompson et al.: https://doi.org/10.1063/1.3245303
                const auto virialI = forceI * i.getR();
                if (i.isOwned()) {
                    _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
                    _aosThreadDataGlobals[threadnum].virialSum += virialI;
                }
                // for non-newton3 particles j and/or k will be considered in a separate calculation
                if (newton3 and j.isOwned()) {
                    const auto virialJ = forceJ * j.getR();
                    _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
                    _aosThreadDataGlobals[threadnum].virialSum += virialJ;
                }
                if (newton3 and k.isOwned()) {
                    const auto virialK = forceK * k.getR();
                    _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy;
                    _aosThreadDataGlobals[threadnum].virialSum += virialK;
                }
                if constexpr (countFLOPs) {
                    if (newton3) {
                        ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
                    } else {
                        ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
                    }
                }
            }
        }




        //added by me for lut
//        [[nodiscard]] std::array<double, 3> getLUTValues(double distSquaredJK,double distSquaredIJ,  double distSquaredKI) const {
        [[nodiscard]] std::array<double, 3> getLUTValues2(double distSquaredIJ,  double distSquaredKI ,double distSquaredJK ) const {

            // Actual distances
            const double distIJ = std::sqrt(distSquaredIJ);
            const double distJK = std::sqrt(distSquaredJK);
            const double distKI = std::sqrt(distSquaredKI);

            // Numerators of cosine representation (cos_i = (r_ij^2 + r_ik^2 - r_jk^2) / (2 * r_ij * r_ik)
            const double numKI = distSquaredIJ + distSquaredJK - distSquaredKI;
            const double numJK = distSquaredIJ + distSquaredKI - distSquaredJK;
            const double numIJ = distSquaredJK + distSquaredKI - distSquaredIJ;

            const double numerator = numKI * numJK * numIJ;

            const double allDistsSquared = distSquaredIJ * distSquaredJK * distSquaredKI;
            const double allDists = distIJ * distJK * distKI;
            const auto allDistsTripled = allDistsSquared * allDists;


            // Gradient factors of 1. / (rrr)^3
            const double allDistsTriplesGradientIJ = 3. / (allDistsTripled * distSquaredIJ);
            const double allDistsTriplesGradientKI = -3. / (allDistsTripled * distSquaredKI);

            // Product of all cosines multiplied with 3: 3 * cos(a)cos(b)cos(c)
           const auto cosines = (3. / 8.) * numerator / allDistsSquared;
            const double cosinesGradientIJ =
                    (3. / 4.) *
                    ((numerator / distSquaredIJ - numKI * numIJ - numJK * numIJ + numJK * numKI) / allDistsSquared);
            const double cosinesGradientKI =
                    (3. / 4.) * ((-numerator / distSquaredKI + numKI * numIJ - numJK * numIJ + numJK * numKI) /
                                 allDistsSquared);

            // Gradient factors corresponding to the normal ATM term
            const auto fullATMGradientIJ =
                    _nu * ((1. + cosines) * allDistsTriplesGradientIJ + cosinesGradientIJ / allDistsTripled);
            const auto fullATMGradientKI =
                    _nu * ((1. + cosines) * allDistsTriplesGradientKI + cosinesGradientKI / allDistsTripled);

            const auto expTerm = std::exp(-_alpha * (distIJ + distJK + distKI));

            // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
            std::array<double, 6> sumFactors{};
            double sum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
                sum += sumFactors[n];
            }

            // Gradient factor of the sum in ij-direction
            double ijSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha);
            }

            // Gradient factor of the sum in ki-direction
            double kiSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha);
            }


            // Total gradient factors for the exponential term times the cosines term
            const double fullExpGradientIJ = expTerm * (-(1. + cosines) * ijSum / distIJ + cosinesGradientIJ * sum);
            const double fullExpGradientKI = expTerm * ((1. + cosines) * kiSum / distKI + cosinesGradientKI * sum);


            const auto factorForceJDirectionIJ = (fullATMGradientIJ + fullExpGradientIJ);
            const auto factorForceJDirectionKI = (fullATMGradientKI + fullExpGradientKI);



            //for the newton part

            // Calculate all components for jk-direction
            const double allDistsTriplesGradientJK = 3. / (allDistsTripled * distSquaredJK);
            const double cosinesGradientJK =
                    (3. / 4.) * ((numerator / distSquaredJK + numKI * numIJ - numJK * numIJ - numJK * numKI) / allDistsSquared);
            const auto fullATMGradientJK =
                    _nu * ((1. + cosines) * allDistsTriplesGradientJK + cosinesGradientJK / allDistsTripled);

            double jkSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
            }
            const double fullExpGradientJK = expTerm * (-(1. + cosines) * jkSum / distJK + cosinesGradientJK * sum);


            const auto factorForceJDirectionJK = (fullATMGradientJK + fullExpGradientJK);

//            return {factorForceJDirectionJK, factorForceJDirectionIJ, factorForceJDirectionKI};
            return { factorForceJDirectionIJ, factorForceJDirectionKI, factorForceJDirectionJK};

        }



        [[nodiscard]] std::array<double, 3> getLUTValues(double dist1Squared,  double dist2Squared ,double dist3Squared ) const {

          double distIJ = std::sqrt(dist1Squared);
          double distJK = std::sqrt(dist2Squared);
          double distKI = std::sqrt(dist3Squared);


            // Calculate prefactor
          const auto allDist = distIJ * distJK *distKI;
            const double allDistsSquared = dist1Squared * dist2Squared * dist3Squared;
            const double allDistsTo5 = allDistsSquared * allDistsSquared * std::sqrt(allDistsSquared);
//            const auto allDistsTripled = allDistsSquared * allDists;

            const auto allDistCubed = (allDist) * (allDist)*(allDist);

            const double factor =  _nu / allDistCubed;  // C_ATM / (R1R2R3)^3


            // Gradient factors of 1. / (rrr)^3
            const double allDistsTriplesGradientIJ = -3. * (distJK* distKI) / (allDist);
            const double allDistsTriplesGradientKI = -3.  * (distJK* distIJ) / (allDist);
            const double allDistsTriplesGradientJJ = -3.  * (distIJ* distIJ) / (allDist);


            // Dot products of both distance vectors going from one particle
            const double IJDotKI = - 0.5 * (dist1Squared + dist3Squared - dist2Squared);
            const double IJDotJK = - 0.5 * (dist1Squared + dist2Squared - dist3Squared);
            const double JKDotKI = - 0.5 * (dist2Squared + dist3Squared - dist1Squared);

            const double allDotProducts = IJDotKI * IJDotJK * JKDotKI;


            const auto forceIDirIJ = IJDotJK * JKDotKI - dist2Squared * dist3Squared + 5.0 * allDotProducts / dist1Squared;
            const auto forceJDirJK = IJDotKI * JKDotKI - dist1Squared * dist3Squared + 5.0 * allDotProducts / dist2Squared;
            const auto forceKDirKI = IJDotJK * JKDotKI - dist1Squared * dist2Squared + 5.0 * allDotProducts / dist3Squared;
            const auto forceIDirJK = IJDotKI * (IJDotJK - JKDotKI);
            const auto forceJDirKI =  IJDotJK * (JKDotKI - IJDotKI);



            double dist_sum = distIJ + distKI + distJK;
            double dist_prod = distIJ * distKI * distJK;

            // Calculate factors and sum for: \sum_{n=0}^5 A_{2n}(r_ij*r_jk*r_ki)^(2n/3)
            std::array<double, 6> sumFactors{};
            double sum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                sumFactors[n] = _constantsA[n] * std::pow((distIJ * distJK * distKI), 2. * n / 3.);
                sum += sumFactors[n];
            }

            const auto expTerm = std::exp(-1* _alpha * (distIJ + distJK + distKI));

            const auto expandedTerm = expTerm * sum;
            const auto A = factor + expandedTerm;

            auto devAatm_IJ = -3 * (_nu   / (allDistCubed * distIJ)) ;
            auto devAatm_KI = -3 * (_nu   / (allDistCubed * distKI)) ;
            auto devAatm_JK = -3 * (_nu   / (allDistCubed * distJK)) ;


            const double numKI = dist1Squared + dist2Squared - dist3Squared;
            const double numJK = dist1Squared + dist2Squared - dist3Squared;
            const double numIJ = dist3Squared + dist2Squared - dist1Squared;

            const double numerator = numKI * numJK * numIJ;
            auto cosines = (3. / 8.) * numerator / allDistsSquared;

            //Derivatives of expandedTerm
            // Gradient factor of the sum in ij-direction
            double ijSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
              ijSum += sumFactors[n] * (2. * n / (3. * distIJ) - _alpha);
            }
            // Gradient factor of the sum in ki-direction
            double kiSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
              kiSum += sumFactors[n] * (2. * n / (3. * distKI) - _alpha);
            }

            // Gradient factor of the sum in jk-direction
            double jkSum = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
              jkSum += sumFactors[n] * (2. * n / (3. * distJK) - _alpha);
            }


            double ijSum_dev = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                ijSum_dev += sumFactors[n] * (-1*  _alpha - (2. * n / (3. * distIJ) ));
            }
            double kiSum_dev = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                kiSum_dev += sumFactors[n] * ( -1*  _alpha -  (2. * n / (3. * distKI) ));
            }

            double jkSum_dev = 0.0;
            for (auto n = 0; n < sumFactors.size(); n++) {
                jkSum_dev += sumFactors[n] * (-1 * _alpha - (2. * n / (3. * distJK) ));
            }



//            auto devCorr_IJ =  -1  * _nu * expTerm* sum + expTerm *

//            auto devIJ = -1*  _alpha*expTerm * sum + expTerm *ijSum;
            auto devIJ =     devAatm_IJ + expTerm *  ijSum_dev;
//            auto devIJ =  (-(1. + cosines) * ijSum / distIJ);
//            auto devKI = -_alpha*expTerm * sum + expTerm *kiSum;
            auto devKI =     devAatm_KI + expTerm *  kiSum_dev;
//            auto devKI =  ((1. + cosines) * kiSum / distKI);
//            auto devJK = -_alpha*expTerm * sum + expTerm *jkSum;
            auto devJK =     devAatm_JK + expTerm *  jkSum_dev;
//            auto devJK =  (-(1. + cosines) * jkSum / distJK);
            return { devIJ, devKI, devJK};



        }

        /**
         * @copydoc autopas::Functor::getNeededAttr()
         */
        constexpr static auto getNeededAttr() {
            return std::array<typename Particle::AttributeNames, 9>{
                    Particle::AttributeNames::id, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
                    Particle::AttributeNames::posZ, Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                    Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId,
                    Particle::AttributeNames::ownershipState};
        }

        /**
         * @copydoc autopas::Functor::getNeededAttr(std::false_type)
         */
        constexpr static auto getNeededAttr(std::false_type) {
            return std::array<typename Particle::AttributeNames, 6>{
                    Particle::AttributeNames::id, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
                    Particle::AttributeNames::posZ, Particle::AttributeNames::typeId,
                    Particle::AttributeNames::ownershipState};
        }

        /**
         * @copydoc autopas::Functor::getComputedAttr()
         */
        constexpr static auto getComputedAttr() {
            return std::array<typename Particle::AttributeNames, 3>{
                    Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                    Particle::AttributeNames::forceZ};
        }

        /**
         * Reset the global values.
         * Will set the global values to zero to prepare for the next iteration.
         */
        void initTraversal() final {
            _potentialEnergySum = 0.;
            _virialSum = {0., 0., 0.};
            _postProcessed = false;
            for (size_t i = 0; i < _aosThreadDataGlobals.size(); ++i) {
                _aosThreadDataGlobals[i].setZero();
            }
        }

        /**
         * Accumulates global values, e.g. potential energy and virial.
         * @param newton3
         */
        void endTraversal(bool newton3) final {
            using namespace autopas::utils::ArrayMath::literals;

            if (_postProcessed) {
                throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
            }
            if (calculateGlobals) {
                // Accumulate potential energy and virial values.
                for (const auto &data: _aosThreadDataGlobals) {
                    _potentialEnergySum += data.potentialEnergySum;
                    _virialSum += data.virialSum;
                }

                // For each interaction, we added the full contribution for all three particles. Divide by 3 here, so that each
                // contribution is only counted once per triplet.
                _potentialEnergySum /= 3.;

                _postProcessed = true;

                AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
                AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
              logToFile(  std::to_string(_potentialEnergySum), "potentialEnergy_KR_100NN");
              logToFile(  std::to_string(_virialSum[0]) +"," +  std::to_string(_virialSum[1]) +"," +  std::to_string(_virialSum[2]), "virial_KR_100NN");

            }
        }

        /**
         * Get the potential Energy.
         * @return the potential Energy
         */
        double getPotentialEnergy() {
            if (not calculateGlobals) {
                throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Trying to get potential energy even though calculateGlobals is false. If you want this functor to calculate "
                        "global "
                        "values, please specify calculateGlobals to be true.");
            }
            if (not _postProcessed) {
                throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Cannot get potential energy, because endTraversal was not called.");
            }
            return _potentialEnergySum;
        }

        /**
         * Get the virial.
         * @return
         */
        double getVirial() {
            if (not calculateGlobals) {
                throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
                        "values, please specify calculateGlobals to be true.");
            }
            if (not _postProcessed) {
                throw autopas::utils::ExceptionHandler::AutoPasException(
                        "Cannot get virial, because endTraversal was not called.");
            }
            return _virialSum[0] + _virialSum[1] + _virialSum[2];
        }

        /**
         * Gets the number of useful FLOPs.
         *
         * @return number of FLOPs since initTraversal() is called.
         */
        [[nodiscard]] size_t getNumFLOPs() const override {
            if constexpr (countFLOPs) {
                const size_t numDistCallsAcc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
                const size_t numKernelCallsN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
                const size_t numKernelCallsNoN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });
                const size_t numGlobalCalcsN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsN3; });
                const size_t numGlobalCalcsNoN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsNoN3; });

                constexpr size_t numFLOPsPerDistanceCall = 0;
                constexpr size_t numFLOPsPerN3KernelCall = 0;
                constexpr size_t numFLOPsPerNoN3KernelCall = 0;
                constexpr size_t numFLOPsPerN3GlobalCalc = 0;
                constexpr size_t numFLOPsPerNoN3GlobalCalc = 0;

                return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
                       numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall +
                       numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
                       numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
            } else {
                // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
                return std::numeric_limits<size_t>::max();
            }
        }

        [[nodiscard]] double getHitRate() const override {
            if constexpr (countFLOPs) {
                const size_t numDistCallsAcc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
                const size_t numKernelCallsN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
                const size_t numKernelCallsNoN3Acc =
                        std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                        [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

                return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
                       (static_cast<double>(numDistCallsAcc));
            } else {
                // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
                return std::numeric_limits<double>::quiet_NaN();
            }
        }

    private:
        /**
         * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
         */
        class AoSThreadDataGlobals {
        public:
            AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}

            void setZero() {
                virialSum = {0., 0., 0.};
                potentialEnergySum = 0.;
            }

            // variables
            std::array<double, 3> virialSum;
            double potentialEnergySum;

        private:
            // dummy parameter to get the right size (64 bytes)
            double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
        };

        /**
         * This class stores internal data for FLOP counters for each thread. Make sure that this data has proper size, i.e.
         * k*64 Bytes!
         * The FLOP count and HitRate are not counted/calculated directly, but through helper counters (numKernelCallsNoN3,
         * numKernelCallsN3, numDistCalls, numGlobalCalcs) to reduce computational cost in the functors themselves and to
         * improve maintainability (e.g. if the cost of a kernel call changes).
         */
        class AoSThreadDataFLOPs {
        public:
            AoSThreadDataFLOPs() : __remainingTo64{} {}

            /**
             * Set all counters to zero.
             */
            void setZero() {
                numKernelCallsNoN3 = 0;
                numKernelCallsN3 = 0;
                numDistCalls = 0;
                numGlobalCalcsN3 = 0;
                numGlobalCalcsNoN3 = 0;
            }

            /**
             * Number of calls to Lennard-Jones Kernel with newton3 disabled.
             * Used for calculating number of FLOPs and hit rate.
             */
            size_t numKernelCallsNoN3 = 0;

            /**
             * Number of calls to Lennard-Jones Kernel with newton3 enabled.
             * Used for calculating number of FLOPs and hit rate.
             */
            size_t numKernelCallsN3 = 0;

            /**
             * Number of distance calculations.
             * Used for calculating number of FLOPs and hit rate.
             */
            size_t numDistCalls = 0;

            /**
             * Counter for the number of times the globals have been calculated with newton3 enabled.
             */
            size_t numGlobalCalcsN3 = 0;

            /**
             * Counter for the number of times the globals have been calculated without newton3 enabled.
             */
            size_t numGlobalCalcsNoN3 = 0;

        private:
            /**
             * dummy parameter to get the right size (64 bytes)
             */
            double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
        };

        // make sure of the size of AoSThreadDataGlobals
        static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
        static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

        const double _cutoffSquared;

        // Parameters of the extended Axilrod-Teller potential for Krypton in Kelvin (K) and Angström (A)
//        const double _nu = 1.61525e6;   // K*A^9
        const double _nu = 1.61525e-3;   // K·nm^9  // K*A^9
//        const double _alpha = 1.378382;  // A^-1
        const double _alpha = 13.78382;  // A^-1  changed to nm


        // Units: {K, K*A^-2, K*A^-4, K*A^-6, K*A^-8, K*A^-10}
        const std::array<double, 6> _constantsA = {-0.3081304e8, -0.3519442e8, 0.4928052e7, -0.2182411e6, 0.343088e4,
                                                   0.0};

        ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

        LUT3B *_lut = nullptr;

        // sum of the potential energy, only calculated if calculateGlobals is true
        double _potentialEnergySum;

        // sum of the virial, only calculated if calculateGlobals is true
        std::array<double, 3> _virialSum;

        // thread buffer for aos
        std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals;
        std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

        // defines whether or whether not the global values are already preprocessed
        bool _postProcessed;


      //helper for evaluation:
      void logToFile(const std::string& message, std::string filename) {
        std::ofstream outFile(filename + ".txt", std::ios::app); // Open in append mode
        if (outFile.is_open()) {
          outFile << message << std::endl;
        } else {
          std::cerr << "Unable to open file for writing." << std::endl;
        }
      }
    };
}  // namespace mdLib
