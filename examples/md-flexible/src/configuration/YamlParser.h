/**
 * @file YamlParser.h
 * @author N. Fottner, D. Martin
 * @date 15.07.2019, 11.04.2023
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include "MDFlexConfig.h"
#include "autopas/utils/NumberSet.h"
#include "objects/Object.h"

/**
 * Parser for input through YAML files.
 */
namespace MDFlexParser::YamlParser {

/**
 * Custom Exception for the Yaml parser that is thrown if an error occured during parsing.
 */
class YamlParserException : public std::exception {
 private:
  const char *message;

 public:
  /**
   * Constructor for YamlParserException
   * @param msg The message of the exception.
   */
  YamlParserException(const char *msg) : message(msg) {}

  /**
   * @return Returns the message of the exception.
   */
  const char *what() const noexcept override { return message; }
};

/**
 * Parses the Input for the simulation from the Yaml File specified in the configuration
 * @param config configuration where the input is stored.
 * @return false if any errors occurred during parsing.
 */
bool parseYamlFile(MDFlexConfig &config);

/**
 * Parses the scalar value of a key in a Object-Node of a YAML-config.
 * @param it YAML-iterator that points to a key of an Object-Node.
 * @param key The key to parse.
 * @param exp The expected value of the key.
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T>
T parseObjectValue(YAML::iterator &it, const char *key, const char *exp);

/**
 * Parses the sequence value of a key in a Object-Node of a YAML-config
 * @param it YAML-iterator that points to a key of an Object-Node.
 * @param key The key to parse.
 * @param exp The expected value of the key.
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T, size_t S>
std::array<T, S> MDFlexParser::YamlParser::parseObjectValue(YAML::iterator &it, const char *key, const char *exp);

/**
 * Parses a CubeGrid-Object from a CubeGrid-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeGrid-Node.
 * @return Particles from the CubeGrid-Generator.
 */
CubeGrid parseCubeGridObject(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses a CubeUniform-Object from a CubeUniform-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeUniform-Node.
 * @return Particles from the CubeUniform-Generator.
 */
CubeUniform parseCubeUniformObject(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses a CubeGauss-Object from a CubeGauss-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeGauss-Node.
 * @return Particles from the CubeGauss-Generator.
 */
CubeGauss parseCubeGaussObject(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses a Sphere-Object from a Sphere-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a Sphere-Node.
 * @return Particles from the Sphere-Generator.
 */
Sphere parseSphereObject(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses a CubeClosestPacked-Object from a CubeClosestPacked-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeClosestPacked-Node.
 * @return Particles from the CubeClosestPacked-Generator.
 */
CubeClosestPacked parseCubeClosestPacked(MDFlexConfig &config, YAML::iterator &it);

/**
 * Throws an exception if a parameter from a particle-generator is missing, or could not be parsed.
 * @param key The YAML-key that is missing, or caused the error.
 * @param exp The expected value of this key.
 * @return
 */
void throwObjectParseException(const char *key, const char *exp);

/**
 * Parses the velcity-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the velocity-node of a particleobject.
 * @return Initial particle velocity in x, y and z direction as array.
 */
std::array<double, 3> parseVelocity(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particle-type-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-type-node of a particleobject.
 * @return Particletype to be generated.
 */
unsigned long parseParticleType(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particle-epsilon-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-epsilon-node of a particleobject.
 * @return Epsilon for particles to be generated.
 */
double parseEpsilon(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particle-sigma-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-sigma-node of a particleobject.
 * @return Sigma for particles to be generated.
 */
double parseSigma(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particle-mass-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-mass-node of a particleobject.
 * @return Mass for particles to be generated.
 */
double parseMass(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particle-spacing-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-spacing-node of a particleobject.
 * @return Spacing between particles to be generated.
 */
double parseParticleSpacing(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the particles-per-dimension-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particles-per-dimension-node of a particleobject.
 * @return Particles per dimension as array.
 */
std::array<unsigned long, 3> parseParticlesPerDim(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the bottomLeftCorner-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the bottomLeftCorner-node of a particleobject.
 * @return Bottom left corner of particles to be generated.
 */
std::array<double, 3> parseBottomLeftCorner(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the distribution-mean-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the distribution-mean-node of a particleobject.
 * @return Mean value for normal distribution in x, y and z direction.
 */
std::array<double, 3> parseDistrMean(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the distribution-stddeviation-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the distribution-stddeviation-node of a particleobject.
 * @return Standard deviation for normal distribution in x, y and z direction as array.
 */
std::array<double, 3> parseDistrStdDev(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the numberOfParticles-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the numberOfParticles-node of a particleobject.
 * @return Number of particles to be generated
 */
size_t parseNumParticles(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the box-length-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the box-length-node of a particleobject.
 * @return Box size in which particles should be generated in x, y and z direction as array.
 */
std::array<double, 3> parseBoxLength(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the center-parameter of a sphere-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the center-node of a sphereobject.
 * @return Center of the sphere object.
 */
std::array<double, 3> parseCenter(MDFlexConfig &config, YAML::iterator &it);

/**
 * Parses the radius-parameter of a sphere-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the radius-node of a sphereobject.
 * @return Radius of the sphere object.
 */
double parseRadius(MDFlexConfig &config, YAML::iterator &it);

}  // namespace MDFlexParser::YamlParser
