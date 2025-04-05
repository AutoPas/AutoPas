/**
* @file DEMParameters.h
* @author Joon Kim
* @date 27/03/2025
*/

#pragma once

namespace demLib {

/**
 * A class that contains all parameters used in DEM functor.
 *
 * Additionally, this class provides getters and setters for the parameters.
 */
class DEMParameters {
 public:
  /**
   * Constructor of the DEMParameters class.
   * Initializes all DEM-specific parameters.
   */
  DEMParameters(double elasticStiffness,
                double normalViscosity,
                double frictionViscosity,
                double rollingViscosity,
                double torsionViscosity,
                double staticFrictionCoeff,
                double dynamicFrictionCoeff,
                double rollingFrictionCoeff,
                double torsionFrictionCoeff,
                double heatConductivity,
                double heatGenerationFactor)
      : _elasticStiffness(elasticStiffness),
        _normalViscosity(normalViscosity),
        _frictionViscosity(frictionViscosity),
        _rollingViscosity(rollingViscosity),
        _torsionViscosity(torsionViscosity),
        _staticFrictionCoeff(staticFrictionCoeff),
        _dynamicFrictionCoeff(dynamicFrictionCoeff),
        _rollingFrictionCoeff(rollingFrictionCoeff),
        _torsionFrictionCoeff(torsionFrictionCoeff),
        _heatConductivity(heatConductivity),
        _heatGenerationFactor(heatGenerationFactor) {}

  /**
   * Destructor of the DEMParameters class.
   */
  ~DEMParameters() = default;

  // Getters
  double getElasticStiffness() const { return _elasticStiffness; }
  double getNormalViscosity() const { return _normalViscosity; }
  double getFrictionViscosity() const { return _frictionViscosity; }
  double getRollingViscosity() const { return _rollingViscosity; }
  double getTorsionViscosity() const { return _torsionViscosity; }
  double getStaticFrictionCoeff() const { return _staticFrictionCoeff; }
  double getDynamicFrictionCoeff() const { return _dynamicFrictionCoeff; }
  double getRollingFrictionCoeff() const { return _rollingFrictionCoeff; }
  double getTorsionFrictionCoeff() const { return _torsionFrictionCoeff; }
  double getHeatConductivity() const { return _heatConductivity; }
  double getHeatGenerationFactor() const { return _heatGenerationFactor; }

 private:
  /**
   * Following attributes are DEM (Discrete Element Method) specific parameters.
   *
   * For detailed explanation of the parameters see: https://mediatum.ub.tum.de/doc/1773224/1773224.pdf
   *
   */
  const double _elasticStiffness;
  const double _normalViscosity;
  const double _frictionViscosity;
  const double _rollingViscosity;
  const double _torsionViscosity;
  const double _staticFrictionCoeff;
  const double _dynamicFrictionCoeff;
  const double _rollingFrictionCoeff;
  const double _torsionFrictionCoeff;
  const double _heatConductivity;
  const double _heatGenerationFactor;
};


}  // namespace demLib