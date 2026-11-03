/**
 * @file ConfigConversion.cpp
 * @date 03.11.2026
 */
#include "ConfigConversion.h"

#include <pybind11/numpy.h>

namespace py = pybind11;

namespace autopasPython {

YAML::Node toYaml(py::handle object, const std::string &path) {
  if (object.is_none()) {
    return YAML::Node{};
  }
  // bool has to come first because it passes as an integer as well.
  if (py::isinstance<py::bool_>(object)) {
    return YAML::Node{object.cast<bool>()};
  }
  if (py::isinstance<py::str>(object)) {
    return YAML::Node{object.cast<std::string>()};
  }
  if (py::isinstance<py::dict>(object)) {
    YAML::Node node{YAML::NodeType::Map};
    for (const auto &[key, value] : object.cast<py::dict>()) {
      const auto keyString = py::str(key).cast<std::string>();
      node[keyString] = toYaml(value, path + "[" + keyString + "]");
    }
    return node;
  }
  if (py::isinstance<py::array>(object)) {
    return toYaml(object.attr("tolist")(), path);
  }
  if (py::isinstance<py::list>(object) or py::isinstance<py::tuple>(object) or py::isinstance<py::set>(object)) {
    YAML::Node node{YAML::NodeType::Sequence};
    size_t index = 0;
    for (const auto &value : object) {
      node.push_back(toYaml(value, path + "[" + std::to_string(index++) + "]"));
    }
    return node;
  }
  // Covers Python integers as well as the scalars of NumPy.
  if (PyIndex_Check(object.ptr())) {
    const auto value = py::reinterpret_steal<py::object>(PyNumber_Index(object.ptr())).cast<long long>();
    return YAML::Node{value};
  }
  if (py::hasattr(object, "__float__")) {
    return YAML::Node{py::float_(py::reinterpret_borrow<py::object>(object)).cast<double>()};
  }
  throw py::type_error(path + " of type " + py::str(py::type::of(object)).cast<std::string>() +
                       " can not be represented in a configuration.");
}

}  // namespace autopasPython
