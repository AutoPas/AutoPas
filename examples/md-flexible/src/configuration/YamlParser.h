/**
 * @file YamlParser.h
 * @author N. Fottner, D. Martin
 * @date 15.07.2019, 11.04.2023
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include <string>

#include "MDFlexConfig.h"
#include "autopas/utils/NumberSet.h"
#include "objects/Object.h"

/**
 * Parser for input through YAML files.
 */
namespace MDFlexParser::YamlParser {

/**
 * Parses the Input for the simulation from the Yaml File specified in the configuration
 * @param config configuration where the input is stored.
 * @return false if any errors occurred during parsing.
 */
bool parseYamlFile(MDFlexConfig &config);

/**
 * Parses an input where exactly one option is expected and throws an error if a sequence of multiple options is given
 * @param node The current YAML::Node, that should be parsed
 * @param errMsg The error message thrown if multiple options are passed
 * @return String representation of the parsed node
 */
const std::string parseSequenceOneElementExpected(const YAML::Node node, std::string errMsg);

/**
 * Parses the scalar value of a key in a Object-Node of a YAML-config.
 * @param it YAML-iterator that points to a key of an Object-Node.
 * @param key The key to parse.
 * @param exp The expected value of the key.
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T>
const T parseObjectValueSingle(YAML::const_iterator &it, std::string key, std::string exp) {
  T value;
  try {
    value = it->second[key].as<T>();
  } catch (const std::exception &e) {
    std::string msg;
    msg.append("Error parsing ");
    msg.append(key);
    msg.append(". Make sure that key \"");
    msg.append(key);
    msg.append("\" exists and has the expected value: ");
    msg.append(exp);
    throw std::runtime_error(msg);
  }
  return value;
}

/**
 * Parses the sequence value of a key in a Object-Node of a YAML-config
 * @param it YAML-iterator that points to a key of an Object-Node.
 * @param key The key to parse.
 * @param exp The expected value of the key.
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T, size_t S>
const std::array<T, S> parseObjectValueSequence(YAML::const_iterator &it, std::string key, std::string exp) {
  std::array<T, S> value;
  try {
    YAML::Node n = it->second[key];
    for (int i = 0; i < S; i++) {
      value[i] = n[i].as<T>();
    }

  } catch (const std::exception &e) {
    std::string msg;
    msg.append("Error parsing ");
    msg.append(key);
    msg.append(". Make sure that key \"");
    msg.append(key);
    msg.append("\" exists and has the expected value: ");
    msg.append(exp);
    throw std::runtime_error(msg);
  }
  return value;
}

/**
 * Parses a CubeGrid-Object from a CubeGrid-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeGrid-Node.
 * @return Particles from the CubeGrid-Generator.
 */
const CubeGrid parseCubeGridObject(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses a CubeUniform-Object from a CubeUniform-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeUniform-Node.
 * @return Particles from the CubeUniform-Generator.
 */
const CubeUniform parseCubeUniformObject(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses a CubeGauss-Object from a CubeGauss-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeGauss-Node.
 * @return Particles from the CubeGauss-Generator.
 */
const CubeGauss parseCubeGaussObject(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses a Sphere-Object from a Sphere-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a Sphere-Node.
 * @return Particles from the Sphere-Generator.
 */
const Sphere parseSphereObject(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses a CubeClosestPacked-Object from a CubeClosestPacked-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the root of a CubeClosestPacked-Node.
 * @return Particles from the CubeClosestPacked-Generator.
 */
const CubeClosestPacked parseCubeClosestPacked(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the velcity-parameter of a particle-generator.
 * @param it YAML-iterator that points to the velocity-node of a particleobject.
 * @return Initial particle velocity in x, y and z direction as array.
 */
const std::array<double, 3> parseVelocity(YAML::const_iterator &it);

/**
 * Parses the particle-type-parameter of a particle-generator.
 * @param it YAML-iterator that points to the particle-type-node of a particleobject.
 * @return Particletype to be generated.
 */
const unsigned long parseParticleType(YAML::const_iterator &it);

/**
 * Parses the particle-epsilon-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-epsilon-node of a particleobject.
 * @return Epsilon for particles to be generated.
 */
const double parseEpsilon(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the particle-sigma-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-sigma-node of a particleobject.
 * @return Sigma for particles to be generated.
 */
const double parseSigma(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the particle-mass-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-mass-node of a particleobject.
 * @return Mass for particles to be generated.
 */
const double parseMass(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the particle-spacing-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particle-spacing-node of a particleobject.
 * @return Spacing between particles to be generated.
 */
const double parseParticleSpacing(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the particles-per-dimension-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the particles-per-dimension-node of a particleobject.
 * @return Particles per dimension as array.
 */
const std::array<unsigned long, 3> parseParticlesPerDim(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the bottomLeftCorner-parameter of a particle-generator.
 * @param it YAML-iterator that points to the bottomLeftCorner-node of a particleobject.
 * @return Bottom left corner of particles to be generated.
 */
const std::array<double, 3> parseBottomLeftCorner(YAML::const_iterator &it);

/**
 * Parses the distribution-mean-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the distribution-mean-node of a particleobject.
 * @return Mean value for normal distribution in x, y and z direction.
 */
const std::array<double, 3> parseDistrMean(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the distribution-stddeviation-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the distribution-stddeviation-node of a particleobject.
 * @return Standard deviation for normal distribution in x, y and z direction as array.
 */
const std::array<double, 3> parseDistrStdDev(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the numberOfParticles-parameter of a particle-generator.
 * @param it YAML-iterator that points to the numberOfParticles-node of a particleobject.
 * @return Number of particles to be generated
 */
const size_t parseNumParticles(YAML::const_iterator &it);

/**
 * Parses the box-length-parameter of a particle-generator.
 * @param config configuration where the input is stored.
 * @param it YAML-iterator that points to the box-length-node of a particleobject.
 * @return Box size in which particles should be generated in x, y and z direction as array.
 */
const std::array<double, 3> parseBoxLength(const MDFlexConfig &config, YAML::const_iterator &it);

/**
 * Parses the center-parameter of a sphere-generator.
 * @param it YAML-iterator that points to the center-node of a sphereobject.
 * @return Center of the sphere object.
 */
const std::array<double, 3> parseCenter(YAML::const_iterator &it);

/**
 * Parses the radius-parameter of a sphere-generator.
 * @param it YAML-iterator that points to the radius-node of a sphereobject.
 * @return Radius of the sphere object.
 */
const double parseRadius(YAML::const_iterator &it);

}  // namespace MDFlexParser::YamlParser
