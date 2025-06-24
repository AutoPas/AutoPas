#pragma once

#include "autopas/containers/P3M/FFT.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/P3M/P3M_traversal.h"

#include <vector>
#include <complex>
#include <array>

namespace autopas {


template <class Particle_T>
class P3M_container : public LinkedCells<Particle_T> {

    using ParticleType = typename LinkedCells<Particle_T>::ParticleType;
    using GridType = typename std::vector<std::vector<std::vector<double>>>;
    using ComplexGridType = std::vector<std::vector<std::vector<std::complex<double>>>>;

    public:
    P3M_container(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, std::array<int, 3> &N, const double cutoff,
              const double skin, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, unsigned int cao = 3)
              : LinkedCells<Particle_T>(boxMin, boxMax, cutoff, skin, rebuildFrequency, cellSizeFactor, loadEstimator) {

        
        grid_dims = N;
        unsigned int maxGridDim = 0;
        for (int i = 0; i < 3; i++){
            //TODO
            //powerOfTwo(grid_dims[i]);
            if(grid_dims[i] > maxGridDim){
                maxGridDim = grid_dims[i];
            }
        }
        fft = FFT(maxGridDim);

        rs_grid = GridType(grid_dims[0]);
        for(int i = 0; i < grid_dims[0]; i++){
            rs_grid[i] = std::vector<std::vector<double>>(grid_dims[1]);
            for(int k = 0; k < grid_dims[1]; k++){
                rs_grid[i][k] = std::vector<double>(grid_dims[2]);
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
        
    }

    /** can be done once at the start of the simulation/tuning, as it only depends on 
     *  the number of gridpoints in every direction, alpha, grid_dist and the charge_assignment_function
     */
    void computeInfluenceFunction(){
        return;
    }

    /**
     * Checks if a given traversal is allowed for P3M and sets it up for the force interactions.
     * @tparam Traversal Traversal type. E.g. pairwise, triwise
     * @param traversal
     */
    template <typename Traversal>
    void prepareTraversal(Traversal &traversal) {
        auto *traversalP3M = dynamic_cast<P3M_traversal<typename LinkedCells<Particle_T>::ParticleCell> *>(traversal);
        if (traversalP3M) {
            traversalP3M->set_traversal_parameters(cao, grid_dims, this->boxMin, grid_dist, rs_grid, rs_grid_shifted, 
                ks_grid, ks_grid_shifted, optForceInfluence, fft, this->begin(autopas::IteratorBehavior::owned));
        } else {
        autopas::utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the P3M container. TraversalID: {}",
          traversal->getTraversalType());
        }
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

    // real_space charge grid
    GridType rs_grid;
    GridType rs_grid_shifted;
    // transformed grid
    ComplexGridType ks_grid;
    ComplexGridType ks_grid_shifted;


    // array for all ks_grid points
    ComplexGridType optForceInfluence;
    // array for all ks_grid points
    ComplexGridType optEnergyInfluence;
};
}