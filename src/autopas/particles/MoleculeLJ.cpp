/**
 * @file MoleculeLJ.cpp
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */
#include "autopas/particles/MoleculeLJ.h"

namespace autopas {

double MoleculeLJ::EPSILON;
double MoleculeLJ::SIGMA;
double MoleculeLJ::MASS;   //@todo Masse und Old_Force des Particle anders implementieren
std::array<double,3>  MoleculeLJ::OLDF;
}  // namespace autopas
