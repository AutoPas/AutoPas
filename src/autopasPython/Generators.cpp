/**
 * @file Generators.cpp
 * @date 03.11.2026
 */
#include "Generators.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <optional>
#include <stdexcept>

#include "src/configuration/objects/CubeClosestPacked.h"
#include "src/configuration/objects/CubeGauss.h"
#include "src/configuration/objects/CubeGrid.h"
#include "src/configuration/objects/CubeUniform.h"
#include "src/configuration/objects/Sphere.h"

namespace py = pybind11;
using py::literals::operator""_a;

namespace autopasPython {
namespace {

using Vec3 = std::array<double, 3>;
constexpr Vec3 zero{0., 0., 0.};

/**
 * Exactly one of the two ways to specify how dense an object is must be given.
 */
void requireExactlyOne(const std::optional<double> &spacing, const std::optional<double> &density, const char *name) {
  if (spacing.has_value() == density.has_value()) {
    throw std::invalid_argument(std::string(name) + " requires exactly one of particle_spacing and density.");
  }
}

}  // namespace

py::dict particlesToArrays(std::vector<ParticleType> particles) {
  std::ranges::sort(particles, {}, [](const auto &p) { return p.getID(); });
  const auto n = static_cast<py::ssize_t>(particles.size());
  py::array_t<double> positions({n, py::ssize_t{3}}), velocities({n, py::ssize_t{3}}), forces({n, py::ssize_t{3}});
  py::array_t<unsigned long> ids(n), typeIds(n);
  auto r = positions.mutable_unchecked<2>();
  auto v = velocities.mutable_unchecked<2>();
  auto f = forces.mutable_unchecked<2>();
  auto id = ids.mutable_unchecked<1>();
  auto typeId = typeIds.mutable_unchecked<1>();
  for (py::ssize_t i = 0; i < n; ++i) {
    const auto &p = particles[i];
    for (py::ssize_t d = 0; d < 3; ++d) {
      r(i, d) = p.getR()[d];
      v(i, d) = p.getV()[d];
      f(i, d) = p.getF()[d];
    }
    id(i) = p.getID();
    typeId(i) = p.getTypeId();
  }
  return py::dict("positions"_a = positions, "velocities"_a = velocities, "forces"_a = forces, "ids"_a = ids,
                  "type_ids"_a = typeIds);
}

void bindGenerators(py::module_ &module) {
  auto generators =
      module.def_submodule("generators", "Generators creating objects made of particles, as known from md-flexible.");

  py::class_<Object, std::shared_ptr<Object>>(generators, "Generator", "Base class of all particle generators.")
      .def_property_readonly("num_particles", &Object::getParticlesTotal)
      .def_property_readonly("box_min", &Object::getBoxMin)
      .def_property_readonly("box_max", &Object::getBoxMax)
      .def(
          "generate",
          [](const Object &object) {
            std::vector<ParticleType> particles{};
            object.generate(particles);
            return particlesToArrays(std::move(particles));
          },
          "Returns the particles as dictionary of NumPy arrays, ready to be used with AutoPas.add_particles().")
      .def("__repr__", [](const Object &object) { return object.getObjectType() + ":\n" + object.to_string(); });

  py::class_<CubeGrid, Object, std::shared_ptr<CubeGrid>>(generators, "CubeGrid", "Particles on a regular grid.")
      .def(py::init([](const std::array<size_t, 3> &particlesPerDimension, std::optional<double> spacing,
                       std::optional<double> density, const Vec3 &bottomLeftCorner, const Vec3 &velocity,
                       size_t typeId, bool centered) {
             requireExactlyOne(spacing, density, "CubeGrid");
             return spacing ? std::make_shared<CubeGrid>(velocity, typeId, particlesPerDimension, *spacing,
                                                         bottomLeftCorner, centered)
                            : std::make_shared<CubeGrid>(velocity, typeId, particlesPerDimension, bottomLeftCorner,
                                                         *density, centered);
           }),
           py::arg("particles_per_dimension"), py::kw_only(), py::arg("particle_spacing") = py::none(),
           py::arg("density") = py::none(), py::arg("bottom_left_corner") = zero, py::arg("velocity") = zero,
           py::arg("type_id") = 0ul, py::arg("centered") = false);

  py::class_<CubeUniform, Object, std::shared_ptr<CubeUniform>>(generators, "CubeUniform",
                                                                "Uniformly distributed particles.")
      .def(py::init([](const Vec3 &boxLength, std::optional<size_t> numParticles, std::optional<double> density,
                       const Vec3 &bottomLeftCorner, const Vec3 &velocity, size_t typeId) {
             if (numParticles.has_value() == density.has_value()) {
               throw std::invalid_argument("CubeUniform requires exactly one of num_particles and density.");
             }
             return numParticles
                        ? std::make_shared<CubeUniform>(velocity, typeId, *numParticles, boxLength, bottomLeftCorner)
                        : std::make_shared<CubeUniform>(velocity, typeId, boxLength, bottomLeftCorner, *density);
           }),
           py::arg("box_length"), py::kw_only(), py::arg("num_particles") = py::none(),
           py::arg("density") = py::none(), py::arg("bottom_left_corner") = zero, py::arg("velocity") = zero,
           py::arg("type_id") = 0ul);

  py::class_<CubeGauss, Object, std::shared_ptr<CubeGauss>>(generators, "CubeGauss",
                                                            "Normally distributed particles.")
      .def(py::init([](size_t numParticles, const Vec3 &boxLength, const Vec3 &mean, const Vec3 &stddev,
                       const Vec3 &bottomLeftCorner, const Vec3 &velocity, size_t typeId) {
             return std::make_shared<CubeGauss>(velocity, typeId, numParticles, boxLength, mean, stddev,
                                                bottomLeftCorner);
           }),
           py::arg("num_particles"), py::arg("box_length"), py::kw_only(), py::arg("mean"), py::arg("stddev"),
           py::arg("bottom_left_corner") = zero, py::arg("velocity") = zero, py::arg("type_id") = 0ul);

  py::class_<Sphere, Object, std::shared_ptr<Sphere>>(generators, "Sphere", "A sphere of particles on a grid.")
      .def(py::init([](const Vec3 &center, int radius, double spacing, const Vec3 &velocity, size_t typeId) {
             return std::make_shared<Sphere>(velocity, typeId, center, radius, spacing);
           }),
           py::arg("center"), py::arg("radius"), py::kw_only(), py::arg("particle_spacing"),
           py::arg("velocity") = zero, py::arg("type_id") = 0ul);

  py::class_<CubeClosestPacked, Object, std::shared_ptr<CubeClosestPacked>>(generators, "CubeClosestPacked",
                                                                            "Particles in a closest packing.")
      .def(py::init([](const Vec3 &boxLength, std::optional<double> spacing, std::optional<double> density,
                       const std::string &structure, const Vec3 &bottomLeftCorner, const Vec3 &velocity,
                       size_t typeId, bool centered) {
             requireExactlyOne(spacing, density, "CubeClosestPacked");
             if (structure != "fcc" and structure != "hcp") {
               throw std::invalid_argument("The structure has to be either 'fcc' or 'hcp'.");
             }
             const auto lattice = structure == "fcc" ? CubeClosestPacked::FCC : CubeClosestPacked::HCP;
             return spacing ? std::make_shared<CubeClosestPacked>(velocity, typeId, *spacing, boxLength,
                                                                  bottomLeftCorner, lattice, centered)
                            : std::make_shared<CubeClosestPacked>(velocity, typeId, boxLength, bottomLeftCorner,
                                                                  *density, lattice, centered);
           }),
           py::arg("box_length"), py::kw_only(), py::arg("particle_spacing") = py::none(),
           py::arg("density") = py::none(), py::arg("structure") = "fcc", py::arg("bottom_left_corner") = zero,
           py::arg("velocity") = zero, py::arg("type_id") = 0ul, py::arg("centered") = false);
}

}  // namespace autopasPython
