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
 * Creates a string representation of types used in particle objects
 * @return String representation of type
 */
template <typename T>
std::string typeToStr() {
  if (typeid(T) == typeid(double)) {
    return "Double";
  }
  if (typeid(T) == typeid(unsigned int)) {
    return "Unsigned Integer";
  }
  if (typeid(T) == typeid(unsigned long)) {
    return "Unsigned Long";
  }
  if (typeid(T) == typeid(int)) {
    return "Integer";
  }
  if (typeid(T) == typeid(bool)) {
    return "Boolean";
  }
  return "Unknown Type";
}

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
const std::string parseSequenceOneElementExpected(const YAML::Node node, const std::string &errMsg);

/**
 * Creates an error message for all non-object-keys
 * @param mark The Yaml-Mark for printing line- and column number
 * @param key The key that caused the error
 * @param errorMsg Message thrown by an exception in yaml-cpp, in the parser, or in AutoPas
 * @param expected The expected value of the key
 * @param description The parameter description of the key
 * @return String representation of the parsed node
 */
const std::string makeErrorMsg(const YAML::Mark &mark, const std::string &key, const std::string &errorMsg,
                               const std::string &expected, const std::string &description);

/**
 * Parses the scalar value of a key in a complex-type node of a YAML-config.
 * @param node root-YAML-node of a complex-type.
 * @param key The key to parse.
 * @param complexTypeErrors Vector to store all errors during parsing of one complex-type node
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T>
const T parseComplexTypeValueSingle(const YAML::Node node, const std::string &key, std::vector<std::string> &complexTypeErrors) {
  T value;
  try {
    value = node[key].as<T>();
  } catch (const std::exception &e) {
    std::stringstream ss;
    ss << "Error parsing " << key << ". Make sure that key \"" << key
       << "\" exists and has the expected value: " << typeToStr<T>();
    complexTypeErrors.push_back(ss.str());
  }
  return value;
}

/**
 * Parses the sequence value of a key in a complex-type node of a YAML-config
 * @param node root-YAML-node of a complex-type.
 * @param key The key to parse.
 * @param complexTypeErrors Vector to store all errors during parsing of one complex-type node
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T, size_t S>
const std::array<T, S> parseComplexTypeValueSequence(const YAML::Node node, const std::string &key,
                                                std::vector<std::string> &complexTypeErrors) {
  std::array<T, S> value;
  try {
    YAML::Node n = node[key];
    for (int i = 0; i < S; i++) {
      value[i] = n[i].as<T>();
    }

  } catch (const std::exception &e) {
    std::stringstream ss;
    ss << "Error parsing " << key << ". Make sure that key \"" << key << "\" exists and has the expected value: "
       << "YAML-sequence of " << std::to_string(S) << " " << typeToStr<T>() << " values.";
    complexTypeErrors.push_back(ss.str());
  }
  return value;
}

/**
 * Parses the sequence value of a key in a complex-type node of a YAML-config. Variant for an unknown sequence size.
 * @param node root-YAML-node of a complex-type.
 * @param key The key to parse.
 * @param complexTypeErrors Vector to store all errors during parsing of one complex-type node
 * @return Parsed value of key. Throws a runtime_error if key could not be parsed.
 */
template <typename T>
const std::vector<T> parseComplexTypeValueSequence(const YAML::Node node, const std::string &key,
                                                std::vector<std::string> &complexTypeErrors) {
  std::vector<T> value;
  try {
    YAML::Node n = node[key];
    const auto vecLength = n.size();
    value.reserve(vecLength);
    for (int i = 0; i < vecLength; i++) {
      value.emplace_back(n[i].as<T>());
    }

  } catch (const std::exception &e) {
    std::stringstream ss;
    ss << "Error parsing " << key << ". Make sure that key \"" << key << "\" exists and has the expected value: "
       << "YAML-sequence of " << typeToStr<T>() << " values.";
    complexTypeErrors.push_back(ss.str());
  }
  return value;
}

/**
 * Parses a CubeGrid-Object from a CubeGrid-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param node root-YAML-node of an Object.
 * @param objectErrors Vector to store all errors during parsing of one object
 * @return Particles from the CubeGrid-Generator.
 */
const CubeGrid parseCubeGridObject(const MDFlexConfig &config, const YAML::Node node,
                                   std::vector<std::string> &objectErrors);

/**
 * Parses a CubeUniform-Object from a CubeUniform-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param node root-YAML-node of an Object.
 * @param objectErrors Vector to store all errors during parsing of one object
 * @return Particles from the CubeUniform-Generator.
 */
const CubeUniform parseCubeUniformObject(const MDFlexConfig &config, const YAML::Node node,
                                         std::vector<std::string> &objectErrors);

/**
 * Parses a CubeGauss-Object from a CubeGauss-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param node root-YAML-node of an Object.
 * @param objectErrors Vector to store all errors during parsing of one object
 * @return Particles from the CubeGauss-Generator.
 */
const CubeGauss parseCubeGaussObject(const MDFlexConfig &config, const YAML::Node node,
                                     std::vector<std::string> &objectErrors);

/**
 * Parses a Sphere-Object from a Sphere-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param node root-YAML-node of an Object.
 * @param objectErrors Vector to store all errors during parsing of one object
 * @return Particles from the Sphere-Generator.
 */
const Sphere parseSphereObject(const MDFlexConfig &config, const YAML::Node node,
                               std::vector<std::string> &objectErrors);

/**
 * Parses a CubeClosestPacked-Object from a CubeClosestPacked-Yaml-Node.
 * @param config configuration where the input is stored.
 * @param node root-YAML-node of an Object.
 * @param objectErrors Vector to store all errors during parsing of one object
 * @return Particles from the CubeClosestPacked-Generator.
 */
const CubeClosestPacked parseCubeClosestPacked(const MDFlexConfig &config, const YAML::Node node,
                                               std::vector<std::string> &objectErrors);

}  // namespace MDFlexParser::YamlParser
