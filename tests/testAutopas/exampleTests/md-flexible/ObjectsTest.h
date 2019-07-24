#pragma once
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"
#include "../../../../examples/md-flexible/Objects.h"
#include "autopas/utils/ArrayUtils.h"


//testet die BoxMin und BoxMax funktionen für alle Objects
class ObjectsTest :public AutoPasTestBase{
public:

    ObjectsTest()
    :AutoPasTestBase(),
    _CGrid{CubeGrid(particlesPerDim,1.,velocity,center)},
    _CGauss{CubeGauss(numParticles,boxlength,5.,2.,velocity,center)},
    _CUniform{CubeUniform(numParticles,boxlength,velocity,center)},
    _Sphere{(Sphere(center,5,1.,0,velocity))} {}



protected:
std::array<double, 3> velocity={0.,0.,0.};
size_t numParticles=100;
std::array<double,3> center={5.,5.,5.};
    std::array<double, 3> boxlength={10.,10.,10.};
    std::array<size_t, 3> particlesPerDim={10,10,10};
    CubeGrid _CGrid;
    CubeGauss _CGauss;
    CubeUniform _CUniform;
    Sphere _Sphere;

};