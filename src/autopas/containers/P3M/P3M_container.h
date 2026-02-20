#pragma once

#include "autopas/containers/P3M/FFT.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/P3M/P3M_traversal.h"
#include "autopas/containers/P3M/P3M_shortRangeFunctor.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

#include "autopas/utils/Timer.h"

#include <vector>
#include <complex>
#include <array>

namespace autopas {


template <class Particle_T>
class P3M_container : public LinkedCells<Particle_T> {

    using LinkedParticleCell = typename LinkedCells<Particle_T>::ParticleCell;
    using ParticleType = typename LinkedCells<Particle_T>::ParticleType;
    using GridType = typename std::vector<double>;
    using ComplexGridType = std::vector<std::complex<double>>;

    public:
    P3M_container(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, std::array<unsigned int, 3> &N, const double cutoff,
              const double skin, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, unsigned int cao = 3, double alpha = 0.5)
              : LinkedCells<ParticleType>(boxMin, boxMax, cutoff, skin, rebuildFrequency, cellSizeFactor, loadEstimator) {

        setupTimer = utils::Timer();
        setupTimer.start();
        
        grid_dims = N;
        
        fft = FFT(N);
    
        for (int i = 0; i < 3; i++){
            box_lengths[i] = boxMax[i] - boxMin[i];
        }

        unsigned int totalLength = grid_dims[0] * grid_dims[1] * grid_dims[2];

        rs_grid = GridType(totalLength);

        ks_grid = ComplexGridType(totalLength);
        optInfluence = ComplexGridType(totalLength);

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
        this->alpha = alpha;
        
        // changes if the influence function is optimized for msq-error of forces or energy
        bool force = false;
        if(force){
            computeInfluenceForce();
        }else{
            computeInfluenceEnergy();
        }

        /*std::cout << "Influence: " << std::endl;
        for (unsigned int x = 0; x < grid_dims[0]; x++){
            for(unsigned int y = 0; y < grid_dims[1]; y++){
                for(unsigned int z = 0; z < grid_dims[2]; z++){
                    unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                    std::cout << optInfluence[index1d] << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }*/
        
        

        selfForceCoeffs = std::vector<std::vector<double>>(3);
        for (int dim = 0; dim < 3; dim++){
            selfForceCoeffs[dim] = std::vector<double>(2);
        }
        computeSelfForceCoeff();
        /*std::cout << "Self Forces: " << std::endl;
        for (int dim = 0; dim < 3; dim++){
            std::cout << "dim " << dim << " ";
            for(int i = 0; i < 2; i++){
                std::cout << selfForceCoeffs[dim][i] << ", ";
            }
            std::cout << std::endl;
        }*/
       setupTimer.stop();

       fftTimer = utils::Timer();
       shortRangeTimer = utils::Timer();
       chargeAssignmentTimer = utils::Timer();
       forceInterpolationTimer = utils::Timer();
    }

    ~P3M_container(){
        std::cout << "FFT total time: " << fftTimer.getTotalTime() << " ns" << std::endl;
        std::cout << "short Range interactions total time: " << shortRangeTimer.getTotalTime() << " ns" << std::endl;
        std::cout << "P3M setup total time: " << setupTimer.getTotalTime() << " ns" << std::endl;
        std::cout << "charge assignment total time: " << chargeAssignmentTimer.getTotalTime() << " ns" << std::endl;
        std::cout << "force interpolation total time: " << forceInterpolationTimer.getTotalTime() << " ns" << std::endl;
    }

    /** Computes the influence function optimized for msq-error of forces
     * 
     *  can be done once at the start of the simulation/tuning, as it only depends on 
     *  the number of gridpoints in every direction, alpha, grid_dist and the charge_assignment_function
     */
    void computeInfluenceForce(){
        AUTOPAS_OPENMP(parallel)
        {
            std::vector<int> brillouinShiftedX = std::vector<int>(grid_dims[0]);
            std::vector<int> brillouinShiftedY = std::vector<int>(grid_dims[1]);
            std::vector<int> brillouinShiftedZ = std::vector<int>(grid_dims[2]);
            
            computeBrillouinShift(brillouinShiftedX, grid_dims[0]);
            computeBrillouinShift(brillouinShiftedY, grid_dims[1]);
            computeBrillouinShift(brillouinShiftedZ, grid_dims[2]);

            AUTOPAS_OPENMP(for schedule(static))
            for(unsigned int iz=0; iz < grid_dims[2]; iz++){
                for(unsigned int iy = 0; iy < grid_dims[1]; iy++){
                    for(unsigned int ix = 0; ix < grid_dims[0]; ix++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(ix, iy, iz, grid_dims);
                        if(ix%(grid_dims[0]/2) == 0 and iy%(grid_dims[1]/2) == 0 and iz%(grid_dims[2]/2) == 0){
                            optInfluence[index1d] = std::complex<double>(0.0);
                        }else{
                            optInfluence[index1d] = computeForceInfluenceAt(brillouinShiftedX[ix], brillouinShiftedY[iy], brillouinShiftedZ[iz]);
                        }
                    }
                }
            }
        }
        return;
    }

    /**
     * returns an array where index i gives the grid-position (as an index) in the first brillouin zone for 1-dimension, 
     * which is centered on the lower left corner of the box, therefore indecies greater than gridSize/2 wrap arround to be negative 
     */
    void computeBrillouinShift(std::vector<int> &shifts, int gridSize){
        shifts[0] = 0;
        for(int i = 1; i <= gridSize/2; i++){
            shifts[i] = i;
            shifts[gridSize - i] = -i;
        }
    }

    /**
     * Computes (sin(arg)/arg)^exponent
     */
    double powsinc(double arg, double exponent){
        if(arg == 0.0){
            // sin(0)/0 = 1 -> pow(1, exponent) = 1
            return 1.0;
        }else{
            return pow((sin(arg) / arg), exponent);
        }
    }

    /**
     * computes the Force Influence, could be combined with computation of energy influence
     * from https://doi.org/10.1063/1.2932253
    */
    std::complex<double> computeForceInfluenceAt(int brillouinX, int brillouinY, int brillouinZ){
        // the reciprocal positions are index/L * 2*pi
        double posX = 2*M_PI*brillouinX / (double)box_lengths[0];
        double posY = 2*M_PI*brillouinY / (double)box_lengths[1];
        double posZ = 2*M_PI*brillouinZ / (double)box_lengths[2];

        double numerator = 0.0;

        double denom1 = 0.0;
        double denom2 = 0.0;

        // adjusts over how many brillouin zones we sum, more zones means higher accuracy of the influence function, but also much more computations
        // could be adjusted until a satisfying accuracy has been reached
        const int brillounZones = 1;

        for(int bx = -brillounZones; bx <= brillounZones; bx++){
            double posXshift = posX + (2*M_PI*bx * grid_dims[0] / box_lengths[0]);
            double arg = (posXshift * grid_dist[0] / 2);

            double UX2 = powsinc(arg, 2*cao);

            for(int by = -brillounZones; by <= brillounZones; by++){
                double posYshift = posY + (2*M_PI*by * grid_dims[1] / box_lengths[1]);
                double arg = (posYshift * grid_dist[1] / 2);
                double UY2 = powsinc(arg, 2*cao);
                
                for(int bz = -brillounZones; bz <= brillounZones; bz++){
                    double posZshift = posZ + (2*M_PI*bz * grid_dims[2] / box_lengths[2]);
                    double arg = (posZshift * grid_dist[2] /2);

                    double UZ2 = powsinc(arg, 2*cao);

                    double posShift2 = posXshift * posXshift + posYshift * posYshift + posZshift * posZshift;

                    
                    double U2 = UX2*UY2*UZ2;
                    // 4*PI is multiplied later

                    denom1 += posShift2 * U2;
                    denom2 += U2;
                    //double posposshift = posX * posXshift + posY * posYshift + posZ * posZshift;
                    numerator += U2 * exp(- (posShift2 / (4*alpha*alpha)));
                }
            }
        }
        
        //double pos2 = posX*posX + posY*posY + posZ*posZ;
        //double denominator = pow(denominatorG(posX, posY, posZ),2) * pos2;
        double denominator = denom1*denom2;
        // 4 PI^2 comes from the fact that it is missing from the positions
        numerator *= 4 * M_PI;

        return std::complex<double>(numerator/denominator);
    }

    /**
     * computes a denominator for the computation of the optimal energy influence function
     */
    double denominatorG(double x, double y, double z){
        double shiftX = 2*M_PI* grid_dims[0] / box_lengths[0];
        double shiftY = 2*M_PI* grid_dims[1] / box_lengths[1];
        double shiftZ = 2*M_PI* grid_dims[2] / box_lengths[2];

        // adjusts over how many brillouin zones we sum, more zones means higher accuracy of the influence function, but also much more computations
        // could be adjusted until a satisfying accuracy has been reached
        int brillounZones = 2;

        double denom = 0;
        for(int bx = -brillounZones; bx <= brillounZones; bx++){
            double arg = ((x + bx*shiftX) * grid_dist[0] / 2);

            double UX2;
            if(arg == 0.0){
                // sin(arg)/arg = 1 -> pow(1, 2*cao) = 1
                UX2 = 1.0;
            }else{
                UX2 = pow((sin(arg) / arg), 2.0*cao);
            }
            for(int by = -brillounZones; by <= brillounZones; by++){
                double arg = ((y + by*shiftY) * grid_dist[1] / 2);
                double UY2;
                if(arg == 0.0){
                    // sin(arg)/arg = 1 -> pow(1, 2*cao) = 1
                    UY2 = 1.0;
                }else{
                    UY2 = pow((sin(arg) / arg), 2.0*cao);
                }
                for(int bz = -brillounZones; bz <= brillounZones; bz++){
                    double arg = ((z + bz*shiftZ) * grid_dist[2] /2);

                    double UZ2;
                    if(arg == 0.0){
                        // sin(arg)/arg = 1 -> pow(1, 2*cao) = 1
                        UZ2 = 1.0;
                    }else{
                        UZ2 = pow((sin(arg) / arg), 2.0*cao);
                    }

                    denom += UX2*UY2*UZ2;
                }
            }
        }

        return denom;
    }

    void computeInfluenceEnergy(){
        AUTOPAS_OPENMP(parallel)
        {
            std::vector<int> brillouinShiftedX = std::vector<int>(grid_dims[0]);
            std::vector<int> brillouinShiftedY = std::vector<int>(grid_dims[1]);
            std::vector<int> brillouinShiftedZ = std::vector<int>(grid_dims[2]);
            
            computeBrillouinShift(brillouinShiftedX, grid_dims[0]);
            computeBrillouinShift(brillouinShiftedY, grid_dims[1]);
            computeBrillouinShift(brillouinShiftedZ, grid_dims[2]);

            AUTOPAS_OPENMP(for schedule(static))
            for(unsigned int iz=0; iz < grid_dims[2]; iz++){
                for(unsigned int iy = 0; iy < grid_dims[1]; iy++){
                    for(unsigned int ix = 0; ix < grid_dims[0]; ix++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(ix, iy, iz, grid_dims);
                        if(ix == 0 and iy == 0 and iz == 0){
                            optInfluence[index1d] = std::complex<double>(0.0);
                        }else{
                            optInfluence[index1d] = computeEnergyInfluenceAt(brillouinShiftedX[ix], brillouinShiftedY[iy], brillouinShiftedZ[iz]);
                        }
                    }
                }
            }
        }
        return;
    }

    /**
     * computes the Energy Influence, could be combined with computation of force influence
     * from https://doi.org/10.1063/1.2816570
    */
    std::complex<double> computeEnergyInfluenceAt(int brillouinX, int brillouinY, int brillouinZ){
        // the reciprocal positions are index/L * 2*pi
        double posX = 2*M_PI*brillouinX / (double)box_lengths[0];
        double posY = 2*M_PI*brillouinY / (double)box_lengths[1];
        double posZ = 2*M_PI*brillouinZ / (double)box_lengths[2];

        double numerator = 0.0;

        // adjusts over how many brillouin zones we sum, more zones means higher accuracy of the influence function, but also much more computations
        // could be adjusted until a satisfying accuracy has been reached
        const int brillounZones = 0;

        double pi = M_PI;

        // only checked to see if works with brillounZones = 0 -> 1 Zone
        for(int bx = -brillounZones; bx <= brillounZones; bx++){
            double posXshift = posX + (2*M_PI*bx * grid_dims[0] / box_lengths[0]);
            double arg = (posXshift * grid_dist[0] / 2);

            double UX2 = powsinc(arg, 2*cao);

            for(int by = -brillounZones; by <= brillounZones; by++){
                double posYshift = posY + (2*M_PI*by * grid_dims[1] / box_lengths[1]);
                double arg = (posYshift * grid_dist[1] / 2);
                double UY2 = powsinc(arg, 2*cao);
                for(int bz = -brillounZones; bz <= brillounZones; bz++){
                    if(abs(bx) + abs(by) + abs(bz) > brillounZones){
                        continue;
                    }

                    double posZshift = posZ + (2*M_PI*bz * grid_dims[2] / box_lengths[2]);
                    double arg = (posZshift * grid_dist[2] /2);

                    double UZ2 = powsinc(arg, 2*cao);

                    double posShift2 = posXshift * posXshift + posYshift * posYshift + posZshift * posZshift;

                    
                    double U2 = UX2*UY2*UZ2;
                    // 4*PI is multiplied later

                    numerator +=  U2 * exp(- (posShift2 / (4*alpha*alpha))) / posShift2;
                }
            }
        }
        
        double denominator = pow(denominatorG(posX, posY, posZ),2);
        // 4 PI^2 comes from the fact that it is missing from the positions
        numerator *= 4 * pi;

        return std::complex<double>(numerator/denominator);
    }

    /**
     * potential calculation used in self force computation
     */
    double Usum(int x, int y, int z, int shiftXextra, int shiftYextra, int shiftZextra){
        double shiftX = 2*M_PI* grid_dims[0] / box_lengths[0];
        double shiftY = 2*M_PI* grid_dims[1] / box_lengths[1];
        double shiftZ = 2*M_PI* grid_dims[2] / box_lengths[2];

        double posX = 2*M_PI*x / (double)box_lengths[0];
        double posY = 2*M_PI*y / (double)box_lengths[1];
        double posZ = 2*M_PI*z / (double)box_lengths[2];

        int brillounZones = 2;

        double sum = 0;
        for(int bx = -brillounZones; bx <= brillounZones; bx++){
            double arg1 = ((posX + bx*shiftX) * grid_dist[0] / 2);
            double arg2 = ((posX + (bx + shiftXextra)*shiftX) * grid_dist[0] / 2);

            double UX1 = powsinc(arg1, cao);
            double UX2 = UX1 * powsinc(arg2, cao);

            for(int by = -brillounZones; by <= brillounZones; by++){
                double arg1 = ((posY + by*shiftY) * grid_dist[1] / 2);
                double arg2 = ((posY + (by + shiftYextra)*shiftY) * grid_dist[1] / 2);

                double UY1 = powsinc(arg1, cao);
                double UY2 = UY1 * powsinc(arg2, cao);
                for(int bz = -brillounZones; bz <= brillounZones; bz++){
                    double arg1 = ((posZ + bz*shiftZ) * grid_dist[2] /2);
                    double arg2 = ((posZ + (bz + shiftZextra)*shiftZ) * grid_dist[2] /2);

                    double UZ1 = powsinc(arg1, cao);
                    double UZ2 = UZ1 * powsinc(arg2, cao);

                    sum += UX2*UY2*UZ2;
                }
            }
        }

        return sum;
    }

    // uses approximation of self foce described in https://doi.org/10.1016/j.cpc.2011.01.026
    void computeSelfForceCoeff(){
        std::vector<int> brillouinShiftedX = std::vector<int>(grid_dims[0]);
        std::vector<int> brillouinShiftedY = std::vector<int>(grid_dims[1]);
        std::vector<int> brillouinShiftedZ = std::vector<int>(grid_dims[2]);
        
        computeBrillouinShift(brillouinShiftedX, grid_dims[0]);
        computeBrillouinShift(brillouinShiftedY, grid_dims[1]);
        computeBrillouinShift(brillouinShiftedZ, grid_dims[2]);

        double vol = box_lengths[0] * box_lengths[1] * box_lengths[2];

        for (int dim = 0; dim < 3; dim++){
            double factor = 2 * M_PI / (vol * grid_dist[dim]);
            for (int coeff = 1; coeff <= 2; coeff++){
                factor *= coeff;

                double selfForceCoeff = 0.0;
                AUTOPAS_OPENMP(parallel for schedule(static) reduction(+: selfForceCoeff))
                for (unsigned int zind = 0; zind < grid_dims[2]; zind++){
                    for (unsigned int yind = 0; yind < grid_dims[1]; yind++){
                        for (unsigned int xind = 0; xind < grid_dims[0]; xind++){
                    
                        
                            if(xind == 0 && yind == 0 && zind == 0){
                                continue;
                            }
                            double precoeff;
                            switch (dim){
                            case 0:
                                precoeff = Usum(brillouinShiftedX[xind], brillouinShiftedY[yind], brillouinShiftedZ[zind], coeff, 0, 0);
                                break;
                            case 1:
                                precoeff = Usum(brillouinShiftedX[xind], brillouinShiftedY[yind], brillouinShiftedZ[zind], 0, coeff, 0);
                                break;
                            case 2:
                                precoeff = Usum(brillouinShiftedX[xind], brillouinShiftedY[yind], brillouinShiftedZ[zind], 0, 0, coeff);
                                break;
                            }
                            unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(xind, yind, zind, grid_dims);
                            selfForceCoeff += real(optInfluence[index1d]) * precoeff;
                        }
                    }   
                }
                // end of parallel for
                selfForceCoeffs[dim][coeff-1] = selfForceCoeff * factor;
            }
        }
    }

    // helps to make the traversal pointer an lvalue
    template <typename Traversal>
    void prepHelp(Traversal *t){
        this->prepareTraversal(t);
    }

    /**
     * Checks if a given traversal is allowed for P3M and sets it up for the force interactions.
     * @tparam Traversal Traversal type. E.g. pairwise, triwise
     * @param traversal
     */
    template <typename Traversal>
    void prepareTraversalP3M(Traversal &traversal, LCC08Traversal<LinkedParticleCell, P3M_shortRangeFunctor<Particle_T>> *shortRangeTraversal) {
        this->prepareTraversal(shortRangeTraversal);
        auto *traversalP3M = dynamic_cast<P3MTraversalInterface<LinkedParticleCell> *>(traversal);
        if (traversalP3M) {
            traversalP3M->set_p3m_traversal_parameters(cao, grid_dims, grid_dist, this->getBoxMin(), selfForceCoeffs, this,
                shortRangeTraversal);
            traversalP3M->set_Timers(&fftTimer, &shortRangeTimer, &chargeAssignmentTimer, &forceInterpolationTimer);
        } else {
        autopas::utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the P3M container. TraversalID: {}",
          traversal->getTraversalType());
        }
    }

    void computeInteractions(TraversalInterface *traversal) override {
        // hard coding to use the P3M_shortRangeFunctor for the traversal
        P3M_shortRangeFunctor<Particle_T> f(alpha, this->getCutoff());

        auto shortRangeTraversal =
        LCC08Traversal<LinkedParticleCell, P3M_shortRangeFunctor<Particle_T>>(
            this->getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getCutoff(),
            this->getCellBlock().getCellLength(), DataLayoutOption::aos, false);
        prepareTraversalP3M(traversal, &shortRangeTraversal);

        traversal->initTraversal();
        traversal->traverseParticles();
        traversal->endTraversal();

        /*for(auto iter = this->begin(); iter.isValid(); ++iter){
            if(iter->isOwned()){
                std::cout << "Particle " << iter->getQ() << " pos : " << iter->getR()[0] << ", " << iter->getR()[1] << ", " << iter->getR()[2] << " final F: " << iter->getF()[0] << ", " << iter->getF()[1] << ", " << iter->getF()[2] << std::endl;
            }
        }*/
    }

    ContainerOption getContainerType() const override{
        return ContainerOption::p3m;
    }

    
    public:
    autopas::FFT fft;
    // real_space charge grid
    GridType rs_grid;
    //GridType rs_grid_shifted;
    // transformed grid
    ComplexGridType ks_grid;
    //ComplexGridType ks_grid_shifted;

    // array for all ks_grid points
    ComplexGridType optInfluence;


    private:
    // ewald splitting parameter
    double alpha;

    const int maxCao = 5;
    int cao; 

    // spacing between grid_points
    std::array<unsigned int, 3> grid_dims;
    std::array<double, 3> grid_dist;
    std::array<double, 3> box_lengths;

    // coefficients to approximate selfForce
    std::vector<std::vector<double>> selfForceCoeffs;

    utils::Timer fftTimer;
    utils::Timer shortRangeTimer;
    utils::Timer chargeAssignmentTimer;
    utils::Timer forceInterpolationTimer;
    utils::Timer setupTimer;
};
}