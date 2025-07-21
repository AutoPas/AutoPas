#pragma once

#include <complex>
#include <vector>
#include <array>
#include <stdio.h>
#include "FFT.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/containers/linkedCells/traversals/LCTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/P3M/P3M_shortRangeFunctor.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"

namespace autopas {


//template <class ParticleCell, class Functor>
template <class ParticleCell>
class P3M_traversal : public LCTraversalInterface, public TraversalInterface {

    using GridType = typename std::vector<std::vector<std::vector<double>>>;
    using ComplexGridType = std::vector<std::vector<std::vector<std::complex<double>>>>;
    using ParticleType = typename ParticleCell::ParticleType;


    public:

    P3M_traversal(/*Functor *functor,*/ const double interactionLength, DataLayoutOption dataLayout, bool useNewton3) 
        : TraversalInterface(dataLayout, useNewton3){
    }

    //P3M_traversal()
    //    :cao(cao), grid_dims(grid_dims), grid_dist(grid_dist), rs_grid(rs_grid), rs_grid_shifted(rs_shifted),
    //     ks_grid(ks_grid), ks_grid_shifted(ks_shifted), optForceInfluence(optForceInfluence){

    //        fft = FFT();
    //}

    void set_traversal_parameters(unsigned int cao, std::array<unsigned int, 3> grid_dims, const std::array<double, 3> &boxMin, std::array<double, 3> grid_dist, GridType &rs_grid, /*GridType &rs_shifted,*/ 
        ComplexGridType &ks_grid, /*ComplexGridType &ks_shifted,*/ ComplexGridType &optForceInfluence, FFT &fft, std::vector<std::vector<double>> &selfForceCoeffs,
        ContainerIterator<ParticleType, true, false> &&beginIter, LCC08Traversal<ParticleCell, P3M_shortRangeFunctor<ParticleType>> *shortRangeTraversal){
        this->cao = cao;
        this->grid_dims = grid_dims;
        this->grid_dist = grid_dist;
        this->rs_grid = rs_grid;
        //this->rs_grid_shifted = rs_shifted;
        this->ks_grid = ks_grid;
        //this->ks_grid_shifted = ks_shifted;
        this->optForceInfluence = optForceInfluence;
        this->fft = fft;
        this->boxMin = boxMin;
        this->selfForceCoeffs = selfForceCoeffs;
        this->beginIter = ContainerIterator<ParticleType, true, false>(beginIter);
        this->shortRangeTraversal = shortRangeTraversal;
    }

    private:
    autopas::FFT fft;

    // ewald splitting parameter
    double alpha;

    const int maxCao = 5;
    unsigned int cao; 

    // spacing between grid_points
    std::array<unsigned int, 3> grid_dims;
    std::array<double, 3> grid_dist;
    std::array<double, 3> boxMin;

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

    std::vector<std::vector<double>> selfForceCoeffs;

    // Container Iterator
    ContainerIterator<ParticleType, true, false> beginIter;
    LCC08Traversal<ParticleCell, P3M_shortRangeFunctor<ParticleType>> *shortRangeTraversal;
        

    // assigns the charges of particles to cao number of points in the rs_grid
    void assignChargeDensities(){
        for (auto iter = beginIter; iter.isValid(); ++iter){
            std::array<double, 3> pPos = iter->getR();
            std::array<int, 3> closestGridpoint = {0,0,0};
            std::vector<double> caFractionsX = std::vector<double>(cao);
            std::vector<double> caFractionsY = std::vector<double>(cao);
            std::vector<double> caFractionsZ = std::vector<double>(cao);

            getChargeAssignmentFractions(pPos, closestGridpoint, caFractionsX, caFractionsY, caFractionsZ);
            // charge assignment according to computed fractions
            for(int xi = 0; xi < cao; xi++){
                for(int yi = 0; yi < cao; yi++){
                    for(int zi = 0; zi < cao; zi++){
                        int xii = xi - (cao/2);
                        int yii = yi - (cao/2);
                        int zii = zi - (cao/2);
                        // 1/V factor may not be needed
                        rs_grid[(closestGridpoint[0] + xii + grid_dims[0]) % grid_dims[0]][(closestGridpoint[1] + yii + grid_dims[1])%grid_dims[1]][(closestGridpoint[2] + zii + grid_dims[2]) % grid_dims[2]]
                            += caFractionsX[xi] * caFractionsY[yi] * caFractionsZ[zi] * iter->getQ() * (1./ (grid_dist[0]*grid_dist[1]*grid_dist[2]));
                    }
                }
            }
        }

    }

    // uses B-splines for which the parameters are tabulated in https://doi.org/10.1063/1.477414
    void chargeFraction(double x, std::vector<double> &fractions){
        if(x > 0.5 || x < -0.5){
            std::cout << "error in charge Fraction" << std::endl;
        }
        //assert(x <= 0.5);
        //assert(x >= -0.5);
        switch(cao){
            case 1:
                fractions[0] = 1.;
                break;
            case 2:
                fractions[0] = 0.5*(1-2*x);
                fractions[1] = 0.5*(1+2*x);
                break;
            case 3:
                fractions[0] = 0.125*(1-4*x+4*x*x);
                fractions[1] = 0.25*(3-4*x*x);
                fractions[2] = 0.125*(1+4*x+4*x*x);
                break;
            case 4:
                fractions[0] = (1./48) * (1-6*x+12*x*x-8*x*x*x);
                fractions[1] = (1./48) * (23-30*x-12*x*x+24*x*x*x);
                fractions[2] = (1./48) * (23+30*x-12*x*x-24*x*x*x);
                fractions[3] = (1./48) * (1+6*x+12*x*x+8*x*x*x);
                break;
            case 5:
                fractions[0] = (1./384) * (1-8*x+24*x*x-32*x*x*x+16*x*x*x*x);
                fractions[1] = (1./96) * (19-44*x+24*x*x+16*x*x*x-16*x*x*x*x);
                fractions[2] = (1./192) * (115-120*x*x+48*x*x*x*x);
                fractions[3] = (1./96) * (19+44*x+24*x*x-16*x*x*x-16*x*x*x*x);
                fractions[4] = (1./384) * (1+8*x+24*x*x+32*x*x*x+16*x*x*x*x);
            default:
                //error unknown charge assignment order
        }
    }

    // computes the derivatives of the charge Fractions
    void cafDerivative(double x, std::vector<double> &fractions){
        //assert(x <= 0.5);
        //assert(x >= -0.5);
        switch(cao){
            case 1:
                fractions[0] = 1.;
                break;
            case 2:
                fractions[0] = -1.;
                fractions[1] = 1.;
                break;
            case 3:
                fractions[0] = x - 0.5;
                fractions[1] = -2 * x;
                fractions[2] = x + 0.5;
                break;
            case 4:
                fractions[0] = -0.5 *x*x + 0.5 * x - (1./8.);
                fractions[1] = 1.5 *x*x - 0.5 * x - (5./8.);
                fractions[2] = -1.5 *x*x - 0.5 * x + (5./8.);
                fractions[3] = 0.5 *x*x + 0.5 * x + (1./8.);
                break;
            case 5:
                fractions[0] = (1./6.) * x*x*x - 0.25 * x*x + 0.125 * x - (1./48);
                fractions[1] = (-2./3.) * x*x*x + 0.5 *x*x + 0.5 * x - (44./96.);
                fractions[2] = x*x*x - (240./192.) * x;
                fractions[3] = (-2./3.) * x*x*x - 0.5 *x*x + 0.5 * x + (44./96.);
                fractions[4] = (1./6.) * x*x*x + 0.25 * x*x + 0.125 * x + (1./48.);
            default:
                //error unknown charge assignment order
        }
    }
    

    // multiplies the ks_grid with an influence function
    void applyInfluenceFunction(){
        for(unsigned int i = 0; i < grid_dims[0]; i++){
            for(unsigned int j = 0; j < grid_dims[1]; j++){
                for(unsigned int k = 0; k < grid_dims[2]; k++){
                    ks_grid[i][j][k] *= optForceInfluence[i][j][k];
                }
            }
        }
    }

    // assigns parts of forces from the points in the rs_grid back to the particles according to the charge assignment fuction
    void interpolateForces(){
        for (auto iter = beginIter; iter.isValid(); ++iter){
            std::array<double, 3> pPos = iter->getR();
            double charge = iter->getQ();
            std::array<int, 3> closestGridpoint = {0,0,0};
            std::vector<double> caFractionsX = std::vector<double>(cao);
            std::vector<double> caFractionsY = std::vector<double>(cao);
            std::vector<double> caFractionsZ = std::vector<double>(cao);
            std::vector<double> caFractionsDX = std::vector<double>(cao);
            std::vector<double> caFractionsDY = std::vector<double>(cao);
            std::vector<double> caFractionsDZ = std::vector<double>(cao);

            getChargeAssignmentFractions(pPos, closestGridpoint, caFractionsX, caFractionsY, caFractionsZ);
            getCAFDeriveative(pPos, closestGridpoint, caFractionsDX, caFractionsDY, caFractionsDZ);
            // charge assignment according to computed fractions
            std::array<double, 3> totalForce = {0., 0., 0.};
            for(int xi = 0; xi < cao; xi++){
                for(int yi = 0; yi < cao; yi++){
                    for(int zi = 0; zi < cao; zi++){
                        int xii = xi - (cao/2);
                        int yii = yi - (cao/2);
                        int zii = zi - (cao/2);

                        double force = rs_grid[(closestGridpoint[0] + xii + grid_dims[0]) % grid_dims[0]][(closestGridpoint[1] + yii + grid_dims[1])%grid_dims[1]][(closestGridpoint[2] + zii + grid_dims[2]) % grid_dims[2]];
                        totalForce[0] -= caFractionsDX[xi] * caFractionsY[yi] * caFractionsZ[zi] * (1./grid_dist[0]) * force * charge;
                        totalForce[1] -= caFractionsX[xi] * caFractionsDY[yi] * caFractionsZ[zi] * (1./grid_dist[1]) * force * charge;
                        totalForce[2] -= caFractionsX[xi] * caFractionsY[yi] * caFractionsDZ[zi] * (1./grid_dist[2]) * force * charge;
                    }
                }
            }
            std::cout << "Particle " << iter->getQ() << " long Range F: " << totalForce[0] << ", " << totalForce[1] << ", " << totalForce[2] << std::endl;
            subtractSelfForce(*iter, totalForce);
            iter->addF(totalForce);
        }
    }

    void subtractSelfForce(ParticleType &p, std::array<double, 3> &force){
        std::cout << " self Force: ";
        for(int dim = 0; dim < 3; dim++){
            double tmp = (selfForceCoeffs[dim][0] * sin(2* M_PI * p.getR()[dim] / grid_dist[dim]) + selfForceCoeffs[dim][1] * sin(4*M_PI*p.getR()[dim] / grid_dist[dim]));
            force[dim] -= p.getQ() * p.getQ() * tmp;
            std::cout << p.getQ() * p.getQ() * tmp << ", ";
        }
        std::cout << std::endl;
    }

    // computes the index of closest Gridpoint in direction dim and saves it in closestGridpoint
    // also returns the position of the closest Gridpoint in direction dim
    double getClosestGridpoint(std::array<double, 3> &pPos,std::array<int, 3> &closestGridpoint, int dim){
        // TODO make more efficient
        double closestGridPos;
        if(cao % 2 == 1){
            int closestPoint = (int)((pPos[dim] - boxMin[dim] + (grid_dist[dim] / 2)) / grid_dist[dim]);
            closestGridpoint[dim] = closestPoint % grid_dims[dim];
            closestGridPos = closestPoint * grid_dist[dim];
        }else{
            int firstPoint = (int)((pPos[dim] - boxMin[dim]) / grid_dist[dim]);
            int secondPoint = (firstPoint + 1);
            closestGridpoint[dim] = secondPoint % grid_dims[dim];
            closestGridPos = (firstPoint * grid_dims[dim] + secondPoint * grid_dims[dim]) / 2;
        }
        return closestGridPos;
    }

    void getChargeAssignmentFractions(std::array<double, 3> &pPos,std::array<int, 3> &closestGridpoint, std::vector<double> &caFractionsX, 
            std::vector<double> &caFractionsY, std::vector<double> &caFractionsZ){
        for(int i = 0; i < 3; i++){
            double closestGridPos = getClosestGridpoint(pPos, closestGridpoint, i);
            switch(i){
                case 0:
                    // TODO currently not in the interval
                    // charge Fraction uses normalized distances in [-0.5, 0.5]
                    chargeFraction((pPos[i] - closestGridPos) / grid_dist[i], caFractionsX);
                    break;
                case 1:
                    chargeFraction((pPos[i] - closestGridPos) / grid_dist[i], caFractionsY);
                    break;
                case 2: 
                    chargeFraction((pPos[i] - closestGridPos) / grid_dist[i], caFractionsZ);
                    break;
            }
            // caFractions now has Charge-fractions
        }
    }

    void getCAFDeriveative(std::array<double, 3> &pPos,std::array<int, 3> &closestGridpoint, std::vector<double> &caFractionsDX, 
            std::vector<double> &caFractionsDY, std::vector<double> &caFractionsDZ){
        for(int i = 0; i < 3; i++){
            double closestGridPos = getClosestGridpoint(pPos, closestGridpoint, i);
            switch(i){
                case 0:
                    // charge Fraction uses normalized distances in [-0.5, 0.5]
                    cafDerivative((pPos[i] - closestGridPos) / grid_dist[i], caFractionsDX);
                    break;
                case 1:
                    cafDerivative((pPos[i] - closestGridPos) / grid_dist[i], caFractionsDY);
                    break;
                case 2: 
                    cafDerivative((pPos[i] - closestGridPos) / grid_dist[i], caFractionsDZ);
                    break;
            }
        }
    }

    /**
     * assigns charges to the rs_grid (2 different shifted grids)
     * transforms the rs_grid
     * applies an influence function to the transformed grid
     * backtransformes the ks_grid
     * interpolates the forces back to the particles (from 2 different girds)
     * 
     * may not fulfill N3
     */
    void traverseFarParticles(){
        assignChargeDensities();
        std::cout << "Charge Grid:" << std::endl;
        for(unsigned int i = 0; i < grid_dims[2]; i++){
            for(unsigned int j = 0; j < grid_dims[1]; j++){
                for(unsigned int k = 0; k < grid_dims[0]; k++){
                    std::cout << rs_grid[k][j][i] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        fft.forward3D(rs_grid, ks_grid, grid_dims);
        applyInfluenceFunction();
         
        fft.backward3D(ks_grid, rs_grid, grid_dims);
        std::cout << "Force Grid:" << std::endl;
        for(unsigned int i = 0; i < grid_dims[2]; i++){
            for(unsigned int j = 0; j < grid_dims[1]; j++){
                for(unsigned int k = 0; k < grid_dims[0]; k++){
                    std::cout << rs_grid[k][j][i] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        interpolateForces();
    }

    /**
     * assigns charges to the rs_grid
     * transforms the rs_grid
     * applies an influence function to the transformed grid
     * for all 3 dimensions:
     * - does ik-differentiation to get forces out
     * - backtransformes the ks_grid
     * - interpolates the forces back to the particles
     * 
     * requires back-transform of vectorial quantities
     */
    /*void traversFarParticlesv2(){
        assignChargeDensities();
        fft.forward();
        applyInfluenceFunction();
        for(dim : dims){
            ikDifferentiation(dim);
            fft.backward();
            interpolateForces(dim);
        }
             
    }*/

    // a speciffic multiplication in the ks_grid
    //void ikDifferentiation();


    // calls normal traversal with modified potential-functor
    void traverseNearParticles(){
        if(shortRangeTraversal){
            shortRangeTraversal->initTraversal();
            shortRangeTraversal->traverseParticles();
            shortRangeTraversal->endTraversal();
            return;
        } 
    }

    public:

    void traverseParticles(){
        traverseFarParticles();
        traverseNearParticles();
    }

    [[nodiscard]] bool isApplicable() const override {
        //utils::isPairwiseFunctor<Functor>()
        return true;//not this->_useNewton3 and this->_dataLayout != DataLayoutOption::soa;
    }

    [[nodiscard]] TraversalOption getTraversalType() const override {
        return TraversalOption::p3m_p3m;
    }

    void initTraversal() override {
        // zero out the grids
        for(unsigned int x = 0; x < grid_dims[0]; x++){
            for(unsigned int y = 0; y < grid_dims[1]; y++){
                for(unsigned int z = 0; z < grid_dims[2]; z++){
                    rs_grid[x][y][z] = 0;
                    ks_grid[x][y][z] = std::complex<double>(0.0);
                }
            }
        }
    }

    void endTraversal() override {
        
    }
};
}