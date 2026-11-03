/**
 * @file ConfigConversion.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/pybind11.h>
#include <yaml-cpp/yaml.h>

#include <string>

namespace autopasPython {

/**
 * Converts dictionaries, sequences, NumPy arrays, and scalars to the YAML node which md-flexible's parser understands.
 * With this a dictionary is nothing but the Python notation of a YAML file, with one parser behind both.
 * @param object The Python object.
 * @param path Where in the configuration this object is. Used for the error message.
 * @return The YAML node.
 * @throws pybind11::type_error if the object has no representation in YAML.
 */
YAML::Node toYaml(pybind11::handle object, const std::string &path = "config");

}  // namespace autopasPython
