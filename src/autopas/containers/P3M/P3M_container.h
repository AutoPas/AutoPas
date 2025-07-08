#pragma once

#include "autopas/containers/P3M/FFT.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/P3M/P3M_traversal.h"
#include "autopas/containers/P3M/P3M_shortRangeFunctor.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"

#include <vector>
#include <complex>
#include <array>

namespace autopas {


template <class Particle_T>
class P3M_container : public LinkedCells<Particle_T> {

    using LinkedParticleCell = LinkedCells<Particle_T>::ParticleCell;
    using ParticleType = typename LinkedCells<Particle_T>::ParticleType;
    using GridType = typename std::vector<std::vector<std::vector<double>>>;
    using ComplexGridType = std::vector<std::vector<std::vector<std::complex<double>>>>;

    public:
    P3M_container(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, std::array<int, 3> &N, const double cutoff,
              const double skin, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, unsigned int cao = 3)
              : LinkedCells<ParticleType>(boxMin, boxMax, cutoff, skin, rebuildFrequency, cellSizeFactor, loadEstimator) {

        
        grid_dims = N;
        unsigned int maxGridDim = 0;
        for (int i = 0; i < 3; i++){
            //TODO move to FFT
            //powerOfTwo(grid_dims[i]);
            if(grid_dims[i] > maxGridDim){
                maxGridDim = grid_dims[i];
            }
        }
        fft = FFT(maxGridDim);
    
        for (int i = 0; i < 3; i++){
            box_lengths[i] = boxMax[i] - boxMin[i];
        } 

        rs_grid = GridType(grid_dims[0]);
        for(int i = 0; i < grid_dims[0]; i++){
            rs_grid[i] = std::vector<std::vector<double>>(grid_dims[1]);
            for(int k = 0; k < grid_dims[1]; k++){
                rs_grid[i][k] = std::vector<double>(grid_dims[2]);
            }
        }

        ks_grid = ComplexGridType(grid_dims[0]);
        optForceInfluence = ComplexGridType(grid_dims[0]);
        for(int i = 0; i < grid_dims[0]; i++){
            ks_grid[i] = std::vector<std::vector<std::complex<double>>>(grid_dims[1]);
            optForceInfluence[i] = std::vector<std::vector<std::complex<double>>>(grid_dims[1]);
            for(int k = 0; k < grid_dims[1]; k++){
                ks_grid[i][k] = std::vector<std::complex<double>>(grid_dims[2]);
                optForceInfluence[i][k] = std::vector<std::complex<double>>(grid_dims[2]);
            }
        }

        for(int i = 0; i < 3; i++){
            grid_dist[i] = (boxMax[i] - boxMin[i]) / N[i];
        }

        if(cao < 1 || cao > maxCao){
            autopas::utils::ExceptionHandler::exception(
             "Unknown Charge assignment of order {}",
             cao);
            return;
        }
        this->cao = cao;
        this->alpha = 0.5;
        
        computeInfluenceFunction();
    }

    /** can be done once at the start of the simulation/tuning, as it only depends on 
     *  the number of gridpoints in every direction, alpha, grid_dist and the charge_assignment_function
     */
    void computeInfluenceFunction(){
        std::vector<int> brillouinShiftedX = std::vector<int>(grid_dims[0]);
        std::vector<int> brillouinShiftedY = std::vector<int>(grid_dims[1]);
        std::vector<int> brillouinShiftedZ = std::vector<int>(grid_dims[2]);
        
        computeBrillouinShift(brillouinShiftedX, grid_dims[0]);
        computeBrillouinShift(brillouinShiftedY, grid_dims[1]);
        computeBrillouinShift(brillouinShiftedZ, grid_dims[2]);

        for(int ix=0; ix < grid_dims[0]; ix++){
            for(int iy = 0; iy < grid_dims[1]; iy++){
                for(int iz = 0; iz < grid_dims[2]; iz++){
                    if(ix%(grid_dims[0]/2) == 0 and iy%(grid_dims[1]/2) == 0 and iz%(grid_dims[2]/2) == 0){
                        optForceInfluence[ix][iy][iz] = std::complex<double>(0.0);
                    }else{
                        optForceInfluence[ix][iy][iz] = computeForceInfluenceAt(brillouinShiftedX[ix], brillouinShiftedY[iy], brillouinShiftedZ[iz]);
                    }
                }
            }
        }
        return;
    }

    void computeBrillouinShift(std::vector<int> shifts, int gridSize){
        for(int i = 0; i < gridSize/2; i++){
            shifts[i] = i;
            shifts[gridSize - i] = -i;
        }
    }

    /*computes the Force Influence, should be combined with computation of energy influence if added*/
    std::complex<double> computeForceInfluenceAt(int brillouinX, int brillouinY, int brillouinZ){
        // the reciprocal positions are index/L * 2*pi
        double posX = brillouinX / (double)box_lengths[0];
        double posY = brillouinY / (double)box_lengths[1];
        double posZ = brillouinZ / (double)box_lengths[2];

        double numerator = 0.0;
        double denominator1 = 0.0;
        double denominator2 = 0.0;

        const int brillounZones = 0;

        for(int bx = -brillounZones; bx <= brillounZones; bx++){
            double posX2 = pow(posX + (bx * grid_dims[0] / box_lengths[0]), 2);
            double arg = (posX * grid_dist[0] + bx) * M_PI;
            double UX2 = pow((sin(arg) / arg), 2.0*cao);
            for(int by = -brillounZones; by <= brillounZones; by++){
                double posY2 = pow(posY + (by * grid_dims[1] / box_lengths[1]), 2);
                double arg = (posY * grid_dist[1] + by) * M_PI;
                double UY2 = pow((sin(arg) / arg), 2.0*cao);
                for(int bz = -brillounZones; bz <= brillounZones; bz++){
                    double posZ2 = pow(posZ + (bz * grid_dims[2] / box_lengths[2]), 2);
                    double arg = (posZ * grid_dist[2] + bx) * M_PI;
                    double UZ2 = pow((sin(arg) / arg), 2.0*cao);

                    double pos2 = posX2 + posY2 + posZ2;

                    double U2 = UX2*UY2*UZ2;

                    // 4 in exp cancels and 4*PI is multiplied later
                    denominator1 += pos2 * U2;
                    denominator2 += U2;
                    numerator += U2 * exp(- (M_PI*M_PI * pos2 / (alpha*alpha)));
                }
            }
        }

        // 4 PI^2 comes from the fact that it is missing from the positions
        double denominator = 4 * M_PI * M_PI * denominator1 * denominator2;
        numerator *= 4 * M_PI;

        return std::complex<double>(numerator/denominator);
    }

    /**
     * Checks if a given traversal is allowed for P3M and sets it up for the force interactions.
     * @tparam Traversal Traversal type. E.g. pairwise, triwise
     * @param traversal
     */
    template <typename Traversal>
    void prepareTraversalP3M(Traversal &traversal, LCC08Traversal<LinkedParticleCell, P3M_shortRangeFunctor<Particle_T>> *shortRangeTraversal) {
        this->prepareTraversal(shortRangeTraversal);
        auto *traversalP3M = dynamic_cast<P3M_traversal<LinkedParticleCell> *>(traversal);
        if (traversalP3M) {
            traversalP3M->set_traversal_parameters(cao, grid_dims, this->getBoxMin(), grid_dist, rs_grid, /*rs_grid_shifted,*/ 
                ks_grid, /*ks_grid_shifted,*/ optForceInfluence, fft, std::move(this->begin(autopas::IteratorBehavior::owned)),
                shortRangeTraversal);
        } else {
        autopas::utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the P3M container. TraversalID: {}",
          traversal->getTraversalType());
        }
    }

    void computeInteractions(TraversalInterface *traversal) override {
        P3M_shortRangeFunctor<Particle_T> f(alpha, this->getCutoff());

        auto shortRangeTraversal =
        LCC08Traversal<LinkedParticleCell, P3M_shortRangeFunctor<Particle_T>>(
            this->getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getCutoff(),
            this->getCellBlock().getCellLength(), DataLayoutOption::aos, false);
        prepareTraversalP3M(traversal, &shortRangeTraversal);

        traversal->initTraversal();
        traversal->traverseParticles();
        traversal->endTraversal();
    }

    ContainerOption getContainerType() const override{
        return ContainerOption::p3m;
    }

    private:
    autopas::FFT fft;

    // ewald splitting parameter
    double alpha;

    const int maxCao = 5;
    int cao; 

    // spacing between grid_points
    std::array<int, 3> grid_dims;
    std::array<double, 3> grid_dist;
    std::array<double, 3> box_lengths;

    // real_space charge grid
    GridType rs_grid;
    //GridType rs_grid_shifted;
    // transformed grid
    ComplexGridType ks_grid;
    //ComplexGridType ks_grid_shifted;


    // array for all ks_grid points
    ComplexGridType optForceInfluence;
    // array for all ks_grid points
    ComplexGridType optEnergyInfluence;
};
}