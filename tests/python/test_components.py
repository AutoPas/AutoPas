"""Tests of the parts which are usable without a simulation, and of the plots."""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from conftest import small_config

import autopas as ap
from autopas import generators
from autopas.plotting import animate_particles, plot_history, plot_particles


@pytest.mark.parametrize(
    "generator, count",
    [
        (generators.CubeGrid([2, 3, 4], particle_spacing=1.0), 24),
        (generators.CubeGrid([2, 2, 2], density=1.0, velocity=[1.0, 0.0, 0.0], type_id=1), 8),
        (generators.CubeUniform([5.0, 5.0, 5.0], num_particles=17), 17),
        (generators.CubeGauss(11, [5.0, 5.0, 5.0], mean=[2.5] * 3, stddev=[1.0] * 3), 11),
    ],
)
def test_generators(generator, count):
    assert generator.num_particles == count
    data = generator.generate()
    assert data["positions"].shape == (count, 3)
    assert data["velocities"].shape == (count, 3)
    np.testing.assert_array_equal(data["ids"], np.arange(count))


def test_sphere_and_closest_packed_are_within_their_bounds():
    for generator in (
        generators.Sphere([5.0, 5.0, 5.0], 3, particle_spacing=1.0),
        generators.CubeClosestPacked([6.0, 6.0, 6.0], particle_spacing=1.1, structure="hcp"),
    ):
        data = generator.generate()
        assert len(data["ids"]) == generator.num_particles > 0
        assert (data["positions"] >= np.array(generator.box_min) - 1e-9).all()
        assert (data["positions"] <= np.array(generator.box_max) + 1e-9).all()


def test_generators_require_one_density_specification():
    with pytest.raises(ValueError, match="exactly one"):
        generators.CubeGrid([2, 2, 2])
    with pytest.raises(ValueError, match="exactly one"):
        generators.CubeGrid([2, 2, 2], particle_spacing=1.0, density=1.0)
    with pytest.raises(ValueError, match="fcc"):
        generators.CubeClosestPacked([2.0] * 3, particle_spacing=1.0, structure="bcc")


def test_generators_in_config_and_as_argument(config):
    first = generators.CubeGrid([2, 2, 2], particle_spacing=1.2, bottom_left_corner=[1.0, 1.0, 1.0])
    second = generators.CubeGrid([2, 2, 2], particle_spacing=1.2, bottom_left_corner=[4.0, 4.0, 4.0], type_id=1)
    sites = {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}, 1: {"mass": 2.0, "epsilon": 1.0, "sigma": 1.0}}
    sim = ap.AutoPas({**config, "Sites": sites, "Objects": [first]}, functors=ap.LennardJones(), objects=[second])
    assert sim.num_particles == 16
    np.testing.assert_array_equal(sim.ids, np.arange(16))
    np.testing.assert_array_equal(np.bincount(sim.type_ids), [8, 8])
    sim(2)


def test_particle_properties_library_mixing():
    library = ap.ParticlePropertiesLibrary(2.5)
    library.add_site_type(0, 1.0)
    library.add_lj_parameters_to_site(0, 1.0, 1.0)
    library.add_site_type(1, 2.0)
    library.add_lj_parameters_to_site(1, 4.0, 2.0)
    library.calculate_mixing_coefficients()
    assert library.num_site_types == 2
    assert library.mass(1) == 2.0
    # Lorentz-Berthelot: epsilon is the geometric mean, sigma the arithmetic one.
    assert library.mixing_24_epsilon(0, 1) == pytest.approx(24.0 * 2.0)
    assert library.mixing_sigma_squared(0, 1) == pytest.approx(1.5**2)


def test_particle_properties_library_in_simulation(config):
    library = ap.ParticlePropertiesLibrary(config["cutoff"])
    library.add_site_type(0, 1.0)
    library.add_lj_parameters_to_site(0, 1.0, 1.0)
    without_sites = {k: v for k, v in config.items() if k != "Sites"}
    sim = ap.AutoPas(without_sites, functors=ap.LennardJones(), particle_properties=library)()
    assert sim.particle_properties.mass(0) == 1.0
    reference = ap.AutoPas(config, functors=ap.LennardJones())()
    np.testing.assert_allclose(sim.forces, reference.forces, atol=1e-9)

    wrong = ap.ParticlePropertiesLibrary(1.0)
    wrong.add_site_type(0, 1.0)
    wrong.add_lj_parameters_to_site(0, 1.0, 1.0)
    with pytest.raises(ValueError, match="cutoff"):
        ap.AutoPas(without_sites, functors=ap.LennardJones(), particle_properties=wrong)


def test_config_accepts_numpy_values(config):
    numpy_config = {
        **config,
        "cutoff": np.float64(2.5),
        "iterations": np.int64(3),
        "box-max": np.array([7.2, 7.2, 7.2]),
    }
    assert ap.AutoPas(numpy_config, functors=ap.LennardJones())().iteration == 3


def test_particle_directly():
    particle = ap.Particle([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 7, 1)
    assert particle.id == 7
    assert particle.type_id == 1
    particle.add_f(np.array([1.0, 1.0, 1.0]))
    particle.sub_f([0.5, 0.5, 0.5])
    np.testing.assert_allclose(particle.f, [0.5, 0.5, 0.5])


def test_plots(config):
    sim = ap.AutoPas(config, functors=ap.LennardJones())(5)
    assert plot_particles(sim).name == "3d"
    assert plot_particles(sim, color_by="speed").name == "3d"
    assert len(plot_history(sim)) == 2
    with pytest.raises(ValueError, match="color_by"):
        plot_particles(sim, color_by="mass")
    plt.close("all")


def test_plot_history_requires_data(config):
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    with pytest.raises(ValueError, match="nothing to plot"):
        plot_history(sim)


@pytest.mark.parametrize("color_by", ["type", "speed"])
def test_animation_advances_the_simulation(config, tmp_path, color_by):
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    animation = animate_particles(sim, frames=3, iterations_per_frame=2, color_by=color_by)
    path = tmp_path / "animation.gif"
    animation.save(path, writer="pillow")
    assert path.stat().st_size > 0
    assert sim.iteration == 6
    plt.close("all")


@pytest.mark.parametrize("speed_max", [None, 5.0])
def test_animation_color_scale_starting_at_rest(config, tmp_path, speed_max):
    config["Objects"]["CubeGrid"][0]["particle-spacing"] = 1.1  # Not in equilibrium, so that particles start to move
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    assert not sim.velocities.any()
    animation = animate_particles(sim, frames=2, iterations_per_frame=2, color_by="speed", speed_max=speed_max)
    animation.save(tmp_path / "animation.gif", writer="pillow")
    scatter = animation._fig.axes[0].collections[0]
    fastest = np.linalg.norm(sim.velocities, axis=1).max()
    assert fastest > 0.0
    assert scatter.get_clim() == (0.0, pytest.approx(fastest if speed_max is None else speed_max))
    plt.close("all")


