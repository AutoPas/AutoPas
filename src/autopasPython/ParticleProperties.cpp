/**
 * @file ParticleProperties.cpp
 * @date 03.11.2026
 */
#include "ParticleProperties.h"

#include <pybind11/stl.h>

#include "src/TypeDefinitions.h"

namespace py = pybind11;

namespace autopasPython {

void bindParticleProperties(py::module_ &module) {
  using PPL = ParticlePropertiesLibraryType;
  py::class_<PPL, std::shared_ptr<PPL>>(module, "ParticlePropertiesLibrary",
                                        "The physical properties of all particle types and their mixed values.")
      .def(py::init<double>(), py::arg("cutoff"))
      .def_property_readonly("cutoff", &PPL::getCutoff)
      .def_property_readonly("num_site_types", &PPL::getNumberRegisteredSiteTypes)
      .def("add_site_type", &PPL::addSiteType, py::arg("site_id"), py::arg("mass"),
           "Registers a type. Types have to be registered with consecutive IDs.")
      .def("add_lj_parameters_to_site", &PPL::addLJParametersToSite, py::arg("site_id"), py::arg("epsilon"),
           py::arg("sigma"))
      .def("add_atm_parameters_to_site", &PPL::addATMParametersToSite, py::arg("site_id"), py::arg("nu"))
      .def("calculate_mixing_coefficients", &PPL::calculateMixingCoefficients,
           "Has to be called after the types changed. The AutoPas constructor does this as well.")
      .def("mass", &PPL::getSiteMass, py::arg("site_id"))
      .def("epsilon", &PPL::getEpsilon, py::arg("site_id"))
      .def("sigma", &PPL::getSigma, py::arg("site_id"))
      .def("nu", &PPL::getNu, py::arg("site_id"))
      .def("mixing_24_epsilon", &PPL::getMixing24Epsilon, py::arg("i"), py::arg("j"))
      .def("mixing_sigma_squared", &PPL::getMixingSigmaSquared, py::arg("i"), py::arg("j"))
      .def("mixing_shift_6", &PPL::getMixingShift6, py::arg("i"), py::arg("j"))
      .def("mixing_nu", &PPL::getMixingNu, py::arg("i"), py::arg("j"), py::arg("k"));
}

}  // namespace autopasPython
