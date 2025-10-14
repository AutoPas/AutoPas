#pragma once

#include <complex>
#include <vector>
#include <array>
#include <stdio.h>
#include <functional>

#include "FFT.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/containers/linkedCells/traversals/LCTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/P3M/P3M_traveralInterface.h"
#include "autopas/containers/P3M/P3M_shortRangeFunctor.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

#include "autopas/utils/Timer.h"

namespace autopas {

// decleration for pointer use later
template <class Particle_Type>
class P3M_container;


template <class ParticleCell, class Functor>
//template <class ParticleCell>
class P3M_traversal : public LCTraversalInterface, public TraversalInterface, public P3MTraversalInterface<ParticleCell> {

    using GridType = typename std::vector<double>;
    using ComplexGridType = std::vector<std::complex<double>>;
    using ParticleType = typename ParticleCell::ParticleType;


    public:

    P3M_traversal(Functor *functor, const double interactionLength, DataLayoutOption dataLayout, bool useNewton3) 
        : TraversalInterface(dataLayout, useNewton3), functor(functor), interactionLength(interactionLength), potential(0.0){
    }

    //P3M_traversal()
    //    :cao(cao), grid_dims(grid_dims), grid_dist(grid_dist), rs_grid(rs_grid), rs_grid_shifted(rs_shifted),
    //     ks_grid(ks_grid), ks_grid_shifted(ks_shifted), optForceInfluence(optForceInfluence){

    //        fft = FFT();
    //}

    void set_p3m_traversal_parameters(unsigned int cao, std::array<unsigned int, 3> grid_dims, std::array<double, 3> grid_dist, const std::array<double, 3> &boxMin,
         std::vector<std::vector<double>> &selfForceCoeffs, P3M_container<ParticleType> *container, LCC08Traversal<ParticleCell, P3M_shortRangeFunctor<ParticleType>> *shortRangeTraversal) override {
        this->cao = cao;
        this->grid_dims = grid_dims;
        this->grid_dist = grid_dist;
        this->boxMin = boxMin;
        this->selfForceCoeffs = selfForceCoeffs;
        this->container = container;
        this->shortRangeTraversal = shortRangeTraversal;
    }

    void set_Timers(utils::Timer *fftTimer, utils::Timer *shortRangeTimer, utils::Timer *chargeAssignmentTimer, utils::Timer *forceInterpolationTimer) override {
        this->fftTimer = fftTimer;
        this->shortRangeTimer = shortRangeTimer;
        this->chargeAssignmentTimer = chargeAssignmentTimer;
        this->forceInterpolationTimer = forceInterpolationTimer;
    }

    private:
    // ewald splitting parameter
    double alpha;

    const int maxCao = 5;
    unsigned int cao; 

    // spacing between grid_points
    std::array<unsigned int, 3> grid_dims;
    std::array<double, 3> grid_dist;
    std::array<double, 3> boxMin;

    std::vector<std::vector<double>> selfForceCoeffs;

    P3M_container<ParticleType> *container;
    LCC08Traversal<ParticleCell, P3M_shortRangeFunctor<ParticleType>> *shortRangeTraversal;
    Functor *functor;
    double interactionLength;
    double potential;

    utils::Timer *fftTimer;
    utils::Timer *shortRangeTimer;
    utils::Timer *chargeAssignmentTimer;
    utils::Timer *forceInterpolationTimer;

    void lamdaFunc (std::vector<double> &caFractionsX, std::vector<double> &caFractionsY, std::vector<double> &caFractionsZ, std::array<int, 3> &closestGridpoint, ParticleType particle) {
        double gridCellVolumeInv = (1./ (grid_dist[0]*grid_dist[1]*grid_dist[2]));

        std::array<double, 3> pPos = particle.getR();

            getChargeAssignmentFractions(pPos, closestGridpoint, caFractionsX, caFractionsY, caFractionsZ);
            // charge assignment according to computed fractions
            for(int zi = 0; zi < cao; zi++){
                int zii = zi - (cao/2);
                unsigned int zIndex = (closestGridpoint[2] + zii + grid_dims[2]) % grid_dims[2];
                for(int yi = 0; yi < cao; yi++){
                    int yii = yi - (cao/2);
                    unsigned int yIndex = (closestGridpoint[1] + yii + grid_dims[1])%grid_dims[1];
                    for(int xi = 0; xi < cao; xi++){
                        int xii = xi - (cao/2);
                        unsigned int xIndex = (closestGridpoint[0] + xii + grid_dims[0]) % grid_dims[0];
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(xIndex, yIndex, zIndex, grid_dims);
                        

                        double chargeFraction = caFractionsX[xi] * caFractionsY[yi] * caFractionsZ[zi] * particle.getQ() * gridCellVolumeInv;
                        AUTOPAS_OPENMP(atomic)
                        container->rs_grid[index1d] += chargeFraction;
                    }
                }
            }
    }

    // assigns the charges of particles to cao number of points in the rs_grid
    void assignChargeDensities(){
        AUTOPAS_OPENMP(parallel){
        double gridCellVolumeInv = (1./ (grid_dist[0]*grid_dist[1]*grid_dist[2]));

        std::array<int, 3> closestGridpoint = {0,0,0};
        std::vector<double> caFractionsX = std::vector<double>(cao);
        std::vector<double> caFractionsY = std::vector<double>(cao);
        std::vector<double> caFractionsZ = std::vector<double>(cao);

        //using region iterators
        int i = autopas_get_thread_num();
        int n = autopas_get_num_threads();
        double iterLength = (container->getBoxMax()[0] - boxMin[0]) / n;
        std::array<double, 3> lowerCorner = {boxMin[0] + iterLength * i, boxMin[1], boxMin[2]};
        std::array<double, 3> upperCorner;
        if(i != n){
            upperCorner = {boxMin[0] + iterLength * (i+1), container->getBoxMax()[1], container->getBoxMax()[2]};
        }else{
            upperCorner = container->getBoxMax();
        }
        
        
        using std::placeholders::_1;
        
        std::function<void(ParticleType)> lambda = std::bind(&P3M_traversal::lamdaFunc, this, caFractionsX, caFractionsY, caFractionsZ, closestGridpoint, _1);
            container->forEachInRegion(lambda , lowerCorner, upperCorner, autopas::IteratorBehavior::owned);
        }

    }

    // uses B-splines for which the parameters are tabulated in https://doi.org/10.1063/1.477414
    void chargeFraction(double x, std::vector<double> &fractions){
        /*if(x > 0.5 || x < -0.5){
            std::cout << "error in charge Fraction" << std::endl;
        }*/
        // TODO remove assert
        assert(x <= 0.5);
        assert(x >= -0.5);
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
                return;
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
                return;
                //error unknown charge assignment order
        }
    }
    

    // multiplies the ks_grid with an influence function
    void applyInfluenceFunction(){
        //compute potential energy
        /*double inv = 1.0 / (grid_dims[0] * grid_dims[1] * grid_dims[2]);
        potential = 0.0;
        AUTOPAS_OPENMP(parallel for schedule(static) reduction(+: potential))
        for(unsigned int i = 0; i < grid_dims[0] * grid_dims[1] * grid_dims[2]; i++){
            potential += inv * inv * std::norm(container->ks_grid[i]) * std::real(container->optForceInfluence[i]);
        }*/

        // *inv is done automatically in the backtransform of the FFT implementation

        AUTOPAS_OPENMP(parallel for schedule(static))
        for(unsigned int i = 0; i < grid_dims[0] * grid_dims[1] * grid_dims[2]; i++){
            container->ks_grid[i] *= container->optForceInfluence[i];
        }
    }

    // assigns parts of forces from the points in the rs_grid back to the particles according to the charge assignment fuction
    void interpolateForces(){
        AUTOPAS_OPENMP(parallel)
        {
        double gridDistInvX = (1./grid_dist[0]);
        double gridDistInvY = (1./grid_dist[1]);
        double gridDistInvZ = (1./grid_dist[2]);

        std::array<int, 3> closestGridpoint = {0,0,0};
        std::vector<double> caFractionsX = std::vector<double>(cao);
        std::vector<double> caFractionsY = std::vector<double>(cao);
        std::vector<double> caFractionsZ = std::vector<double>(cao);
        std::vector<double> caFractionsDX = std::vector<double>(cao);
        std::vector<double> caFractionsDY = std::vector<double>(cao);
        std::vector<double> caFractionsDZ = std::vector<double>(cao);

        for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter){
            std::array<double, 3> pPos = iter->getR();
            double charge = iter->getQ();
            

            getChargeAssignmentFractions(pPos, closestGridpoint, caFractionsX, caFractionsY, caFractionsZ);
            getCAFDeriveative(pPos, closestGridpoint, caFractionsDX, caFractionsDY, caFractionsDZ);
            // charge assignment according to computed fractions
            std::array<double, 3> totalForce = {0., 0., 0.};
            for(int zi = 0; zi < cao; zi++){
                int zii = zi - (cao/2);
                unsigned int zIndex = (closestGridpoint[2] + zii + grid_dims[2]) % grid_dims[2];
                for(int yi = 0; yi < cao; yi++){
                    int yii = yi - (cao/2);
                    unsigned int yIndex = (closestGridpoint[1] + yii + grid_dims[1])%grid_dims[1];
                    for(int xi = 0; xi < cao; xi++){
                        int xii = xi - (cao/2);
                        unsigned int xIndex = (closestGridpoint[0] + xii + grid_dims[0]) % grid_dims[0];
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(xIndex, yIndex, zIndex, grid_dims);
                        
                        

                        double force = container->rs_grid[index1d];
                        totalForce[0] -= caFractionsDX[xi] * caFractionsY[yi] * caFractionsZ[zi] * gridDistInvX * force * charge;
                        totalForce[1] -= caFractionsX[xi] * caFractionsDY[yi] * caFractionsZ[zi] * gridDistInvY * force * charge;
                        totalForce[2] -= caFractionsX[xi] * caFractionsY[yi] * caFractionsDZ[zi] * gridDistInvZ * force * charge;
                    }
                }
            }
            /*AUTOPAS_OPENMP(critical){
            std::cout << "Particle " << iter->getQ() << " long Range F: " << totalForce[0] << ", " << totalForce[1] << ", " << totalForce[2] << std::endl;
            }*/
           //double scaling = 992.573;
          /*totalForce[0] *= scaling;
            totalForce[1] *= scaling;
            totalForce[2] *= scaling;*/
            subtractSelfForce(*iter, totalForce);
            
            iter->addF(totalForce);
        }
    }
    }

    void subtractSelfForce(ParticleType &p, std::array<double, 3> &force){
        //std::cout << " self Force: ";
        for(int dim = 0; dim < 3; dim++){
            double tmp = (selfForceCoeffs[dim][0] * sin(2* M_PI * p.getR()[dim] / grid_dist[dim]) + selfForceCoeffs[dim][1] * sin(4*M_PI*p.getR()[dim] / grid_dist[dim]));
            force[dim] -= p.getQ() * p.getQ() * tmp;
            //std::cout << p.getQ() * p.getQ() * tmp << ", ";
        }
        //std::cout << std::endl;
    }

    // computes the index of closest Gridpoint in direction dim and saves it in closestGridpoint
    // also returns the position of the closest Gridpoint in direction dim
    double getClosestGridpoint(std::array<double, 3> &pPos,std::array<int, 3> &closestGridpoint, int dim){
        double closestGridPos;
        if(cao % 2 == 1){
            int closestPoint = (int)((pPos[dim] - boxMin[dim]) / grid_dist[dim] + 0.5);
            closestGridpoint[dim] = closestPoint % grid_dims[dim];
            closestGridPos = closestPoint * grid_dist[dim];
        }else{
            int firstPoint = (int)((pPos[dim] - boxMin[dim]) / grid_dist[dim]);
            int secondPoint = (firstPoint + 1);
            closestGridpoint[dim] = secondPoint % grid_dims[dim];
            closestGridPos = (firstPoint * grid_dist[dim] + secondPoint * grid_dist[dim]) / 2;
        }
        return closestGridPos;
    }

    void getChargeAssignmentFractions(std::array<double, 3> &pPos,std::array<int, 3> &closestGridpoint, std::vector<double> &caFractionsX, 
            std::vector<double> &caFractionsY, std::vector<double> &caFractionsZ){
        for(int i = 0; i < 3; i++){
            double closestGridPos = getClosestGridpoint(pPos, closestGridpoint, i);
            switch(i){
                case 0:
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
        chargeAssignmentTimer->start();
        assignChargeDensities();
        chargeAssignmentTimer->stop();
        /*std::cout << "Charge Grid:" << std::endl;
        for(unsigned int i = 0; i < grid_dims[2]; i++){
            for(unsigned int j = 0; j < grid_dims[1]; j++){
                for(unsigned int k = 0; k < grid_dims[0]; k++){
                    unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(k, j, i, grid_dims);
                    std::cout << container->rs_grid[index1d] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }*/
        fftTimer->start();
        container->fft.forward3D(container->rs_grid, container->ks_grid, grid_dims);
        fftTimer->stop();

        /*std::cout << "Transformed Charge Density: " << std::endl;
        for (unsigned int x = 0; x < grid_dims[0]; x++){
            for(unsigned int y = 0; y < grid_dims[1]; y++){
                for(unsigned int z = 0; z < grid_dims[2]; z++){
                    unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                    std::cout << container->ks_grid[index1d] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }*/

        applyInfluenceFunction();
        
        fftTimer->start();
        container->fft.backward3D(container->ks_grid, container->rs_grid, grid_dims);
        fftTimer->stop();
        /*std::cout << "Force Grid:" << std::endl;
        for(unsigned int i = 0; i < grid_dims[2]; i++){
            for(unsigned int j = 0; j < grid_dims[1]; j++){
                for(unsigned int k = 0; k < grid_dims[0]; k++){
                    unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(k, j, i, grid_dims);
                    std::cout << container->rs_grid[index1d] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }*/
        forceInterpolationTimer->start();
        interpolateForces();
        forceInterpolationTimer->stop();
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
        }

        auto LJTraversal =
        LCC08Traversal<ParticleCell, Functor>(
            container->getCellBlock().getCellsPerDimensionWithHalo(), functor, container->getCutoff(),
            container->getCellBlock().getCellLength(), DataLayoutOption::aos, false);

        //dynamic_cast<LCC08Traversal<ParticleCell, Functor>>()

        container->prepHelp(&LJTraversal);
        LJTraversal.initTraversal();
        LJTraversal.traverseParticles();
        LJTraversal.endTraversal();
        // LJ functor einfügen
        //TODO Thermostat.h addBrownianMotion() z-coordiante auf 0 setzen.
    }

    public:

    void traverseParticles() override{
        traverseFarParticles();

        shortRangeTimer->start();
        traverseNearParticles();
        shortRangeTimer->stop();
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
        AUTOPAS_OPENMP(parallel for schedule(static))
        for(unsigned int i = 0; i < grid_dims[0] * grid_dims[1] * grid_dims[2]; i++){
            container->rs_grid[i] = 0;
            container->ks_grid[i] = std::complex<double>(0.0);
        }
    }

    void endTraversal() override {
        
    }

    double getPotential() override {
        return potential;
    }
};
}