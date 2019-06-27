//
// Created by nicola on 12.05.19.
//
#pragma once
#include <autopas/utils/MemoryProfiler.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "MDFlexParser.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "TimeDiscretization.h"
#include "autopas/AutoPas.h"


using namespace autopas;
using namespace std;

template<class Particle,class ParticleCell>
class  Simulation {

private:
    shared_ptr<AutoPas<Particle, ParticleCell>> _autopas;

  // autopas::LJFunctor<Particle,ParticleCell, autopas::FunctorN3Modes::Both, true>* _Functor = new
  // autopas::LJFunctor<Particle,ParticleCell, autopas::FunctorN3Modes::Both, true>(1, 1.0, 1.0, 0.0,{0., 0.,
  // 0.},{5., 5., 5.},true);
  long durationX;
  long durationF;
  long durationV;
  double delta_t;
  int iterations;
  ParticleClassLibrary *_PCL;

 public:
  virtual ~Simulation() {}

  explicit Simulation(shared_ptr<AutoPas<Particle, ParticleCell>> autopas, ParticleClassLibrary &PCL);

    /**
    * Writes a VTK file for the current state of the AutoPas object
    * @tparam AutoPasTemplate Template for the templetized autopas type.
    * @param filename
    * @param numParticles
    * @param autopas
    */
    template <class AutoPasTemplate>
    void writeVTKFile(int iteration, size_t numParticles, AutoPasTemplate &autopas) {
        string filename="VtkOutput";
        stringstream strstr;
        strstr << filename << "_" << setfill('0') << setw(4) << iteration << ".vtu";
        //string path = "./vtk";
        std::ofstream vtkFile;
        vtkFile.open(strstr.str());

        vtkFile << "# vtk DataFile Version 2.0" << endl;
        vtkFile << "Timestep" << endl;
        vtkFile << "ASCII" << endl;
        vtkFile << "DATASET STRUCTURED_GRID" << endl;
        vtkFile << "DIMENSIONS 1 1 1" << endl;
        vtkFile << "POINTS " << numParticles << " double" << endl;

        for (auto iter = autopas->begin(); iter.isValid(); ++iter) {
            auto pos = iter->getR();
            vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
        }

        vtkFile.close();
    }

    /**
    * @brief Constructs a container and fills it with particles.
    *
    * According to the options passed, a %DirectSum or %'LinkedCells' container is
    * built. It consists of %`FullParticleCells` and is filled with
    * `PrintableMolecules`. The particles are aligned on a cuboid grid.
    *
    * @param autopas AutoPas object that should be initialized
    * @param particlesPerDim Number of desired particles per dimension.
    * @param particelSpacing Space between two particles along each axis of space.
    */
    void initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas,
                           size_t particlesPerDim, double particelSpacing);

    void initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas,
                            double boxLength, size_t numParticles, double distributionMean, double distributionStdDev);

    void initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas,
                              double boxLength, size_t numParticles);
    /**Getter for Duration of Position Calculation
     * @return durationX
     */
    long getDurationX() const;
    /**Adder for durationX
     * @param durationX
     */
    void addDurationX(long durationX);
    /**Getter for Duration of Force Calculation
     * @return durationF
     */
    long getDurationF() const;
    /**Adder for durationF
     * @param durationF
     */
    void addDurationF(long durationF);
    /**Getter for Duration of Velocity Calculation
     * @return durationV
     */
    long getDurationV() const;
    /**Adder for durationV
     * @param durationV
     */
    void addDurationV(long durationV);

  /** @brief This function
   * -initializes the autopas Object
   * -sets/initializes the simulation domain with the particles generators
   * @todo -initialized Velocities and Positions (and forces?)
   */
  void initialize(MDFlexParser *parser);

    /**
     * Does the ForceCalculation
     * @param Force Calculation Functor
     * @return Duration of Calculation
     * */
    void CalcF();

    /**
     * This function processes the main simulation loop
     * -calls the time discretization class(calculate fores, etc ...)
     * -do the output each timestep
     * -collects the duration of every Calculation(Position,Force,Velocity)
     * @param duration of Simulation
     */
    long simulate();

    /**Getter for Autopas Oject
     * @return Autopas Object
     */
    AutoPas<Particle, ParticleCell> *getAutopas() const;

};
//@TODO add boolean object to choose between LJFunctor and LJF_AVX  -> implement LJF_AVX with mixing rules
template <class Particle, class ParticleCell>
Simulation<Particle, ParticleCell>::Simulation(shared_ptr<AutoPas<Particle, ParticleCell>> autopas,
                                               ParticleClassLibrary &PCL) {
  _autopas = autopas;
  durationF = 0;
  durationV = 0;
  durationX = 0;
  _PCL = &PCL;
}

template<class Particle, class ParticleCell>
AutoPas<Particle, ParticleCell> *Simulation<Particle, ParticleCell>::getAutopas() const {
    return _autopas;
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::initialize(MDFlexParser *parser) {
  // werte die man später für die initialisierung der Funktoren braucht, temporäre implementierung
  // std::array<double, 3> lowCorner = {0., 0., 0.};
  // std::array<double, 3> highCorner = {5., 5., 5.};
  // double epsilon,sigma  = 1.0;
  string logFileName(parser->getLogFileName());
  auto measureFlops(parser->getMeasureFlops());  //@todo un-used
  auto numIterations(parser->getIterations());
  auto particlesTotal(parser->getParticlesTotal());
  auto verletRebuildFrequency(parser->getVerletRebuildFrequency());
  auto vtkFilename(parser->getWriteVTK());  //
  auto boxLength(parser->getBoxLength());
  auto containerChoice(parser->getContainerOptions());
  auto selectorStrategy(parser->getSelectorStrategy());
  auto cutoff(parser->getCutoff());
  auto dataLayoutOptions(parser->getDataLayoutOptions());
  auto distributionMean(parser->getDistributionMean());
  auto distributionStdDev(parser->getDistributionStdDev());
  auto functorChoice(parser->getFunctorOption());  //@todo un-used
  auto generatorChoice(parser->getGeneratorOption());
  auto newton3Options(parser->getNewton3Options());
  auto particleSpacing(parser->getParticleSpacing());
  auto particlesPerDim(parser->getParticlesPerDim());
  auto traversalOptions(parser->getTraversalOptions());
  auto tuningInterval(parser->getTuningInterval());
  auto tuningSamples(parser->getTuningSamples());
  auto verletSkinRadius(parser->getVerletSkinRadius());
  this->delta_t = parser->getDeltaT();
  this->iterations = numIterations;
  _autopas->setCutoff(cutoff);
  _autopas->setVerletSkin(verletSkinRadius);
  _autopas->setVerletRebuildFrequency(verletRebuildFrequency);
  _autopas->setTuningInterval(tuningInterval);
  _autopas->setNumSamples(tuningSamples);
  _autopas->setSelectorStrategy(selectorStrategy);
  _autopas->setAllowedContainers(containerChoice);
  _autopas->setAllowedTraversals(traversalOptions);
  _autopas->setAllowedDataLayouts(dataLayoutOptions);
  _autopas->setAllowedNewton3Options(newton3Options);

    switch (generatorChoice) {
        case MDFlexParser::GeneratorOption::grid: {
            this->initContainerGrid(*_autopas, particlesPerDim, particleSpacing); //particlesTotal wird in diesem fall in der main geupdated
            break;
        }
        case MDFlexParser::GeneratorOption::uniform: {
            this->initContainerUniform(*_autopas, boxLength, particlesTotal);
            break;
        }
        case MDFlexParser::GeneratorOption::gaussian: {
            this->initContainerGauss(*_autopas, boxLength, particlesTotal, distributionMean, distributionStdDev);
            break;
        }
        default:
            throw std::runtime_error("Unknown generator choice");
    }
    this->CalcF();

}
template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGrid(autopas::AutoPas<Particle, ParticleCell> &autopas,
                       size_t particlesPerDim, double particelSpacing) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax(
            {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                     {particelSpacing, particelSpacing, particelSpacing},
                                     {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerGauss(autopas::AutoPas<Particle, ParticleCell> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

template<class Particle,class ParticleCell>
void Simulation<Particle, ParticleCell>::initContainerUniform(autopas::AutoPas<Particle, ParticleCell> &autopas,
                          double boxLength, size_t numParticles) {
    std::array<double, 3> boxMin({0., 0., 0.});
    std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

    autopas.setBoxMin(boxMin);
    autopas.setBoxMax(boxMax);

    autopas.init();

    PrintableMolecule dummyParticle;
    RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

template <class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::CalcF() {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
  //@ TODO: switch for other functors --> mit boolean object?

  // auto* functor = new autopas::LJFunctor<Particle,ParticleCell, autopas::FunctorN3Modes::Both,
  // true>(_autopas->getCutoff(),1., 1.0, 0.0,_autopas->getBoxMin(),_autopas->getBoxMax(),true);
  auto *functor = new LJFunctor<Particle, ParticleCell>(_autopas->getCutoff(), *_PCL, 0.0, _autopas->getBoxMin(),
                                                         _autopas->getBoxMax());
  _autopas->iteratePairwise(functor);
  delete functor;
  functor = NULL;
  stopCalc = std::chrono::high_resolution_clock::now();
  auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  this->addDurationF(durationCalc);
}


template<class Particle,class ParticleCell>
long Simulation<Particle, ParticleCell>::simulate(){
    std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
    startCalc = std::chrono::high_resolution_clock::now();
    double particleDelta_T=this->delta_t;
    double time=0;
    double timeEnd= this->delta_t * (double)this->iterations;
    TimeDiscretization<decltype(_autopas)> td(particleDelta_T);
    //main simulation loop
    while(time<timeEnd){
        this->addDurationX(td.VSCalculateX(_autopas));
        if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
             cout << "Iteration " << time/particleDelta_T << endl;
             cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
        }
        this->CalcF();
        this->addDurationV(td.VSCalculateV(_autopas));
        time+=particleDelta_T;
        this->writeVTKFile<decltype(_autopas)>(time/particleDelta_T,_autopas->getNumberOfParticles(),_autopas);
    }
    stopCalc = std::chrono::high_resolution_clock::now();
    auto durationCalc = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
    return durationCalc;
}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationX() const {
    return durationX;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationX(long durationX) {
    Simulation::durationX += durationX;
}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationF() const {
    return durationF;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationF(long durationF) {
    Simulation::durationF += durationF;
}

template<class Particle, class ParticleCell>
long Simulation<Particle, ParticleCell>::getDurationV() const {
    return durationV;
}

template<class Particle, class ParticleCell>
void Simulation<Particle, ParticleCell>::addDurationV(long durationV) {
    Simulation::durationV += durationV;
}

