"""Tests of constructing and running the AutoPas object."""

import threading
from pathlib import Path

import numpy as np
import pytest
from conftest import CUTOFF, PythonLennardJones, small_config

import autopas as ap

REPO_ROOT = Path(__file__).resolve().parents[2]
NATIVE_FUNCTORS = [ap.LennardJones, ap.LennardJonesHWY]


def assert_same_state(a, b, atol=1e-9):
    np.testing.assert_array_equal(a.ids, b.ids)
    np.testing.assert_allclose(a.positions, b.positions, atol=atol)
    np.testing.assert_allclose(a.velocities, b.velocities, atol=atol)
    np.testing.assert_allclose(a.forces, b.forces, atol=atol)


def test_dict_config_runs_the_configured_iterations(config):
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    assert sim.iteration == 0
    assert sim.num_particles == 216
    assert sim() is sim
    assert sim.iteration == 20
    # Nothing to do anymore, until the configured end is exceeded on request.
    sim()
    assert sim.iteration == 20


def test_resuming_equals_running_at_once(config):
    at_once = ap.AutoPas(config, functors=ap.LennardJones())
    at_once(40)
    in_parts = ap.AutoPas(config, functors=ap.LennardJones())
    in_parts()
    in_parts(5)
    in_parts(15)
    assert in_parts.iteration == 40
    assert_same_state(at_once, in_parts)


def test_yaml_file_equals_dict(tmp_path):
    yaml = """
container: [LinkedCells]
traversal: [lc_c08]
data-layout: [AoS]
newton3: [enabled]
cutoff: 2.5
verlet-skin-radius: 0.3
verlet-rebuild-frequency: 10
deltaT: 0.002
iterations: 20
boundary-type: [periodic, periodic, periodic]
box-min: [0, 0, 0]
box-max: [7.2, 7.2, 7.2]
functor: lj
log-level: warn
Sites:
  0: {mass: 1.0, epsilon: 1.0, sigma: 1.0}
Objects:
  CubeGrid:
    0:
      particles-per-dimension: [6, 6, 6]
      particle-spacing: 1.2
      bottomLeftCorner: [0.1, 0.1, 0.1]
      velocity: [0, 0, 0]
      particle-type-id: 0
"""
    path = tmp_path / "config.yaml"
    path.write_text(yaml)
    from_file = ap.AutoPas(path)
    assert [f.name for f in from_file.functors] == ["LennardJones"]
    from_file()
    from_dict = ap.AutoPas(small_config(), functors=ap.LennardJones())()
    assert_same_state(from_file, from_dict)


@pytest.mark.parametrize("name", ["explodingLiquid.yaml", "fallingDrop.yaml"])
def test_shipped_yaml_files_can_be_loaded(name, tmp_path, monkeypatch):
    # These files ask for VTK output, whose folder is created relative to the working directory.
    monkeypatch.chdir(tmp_path)
    sim = ap.AutoPas(str(REPO_ROOT / "examples/md-flexible/input" / name))
    assert sim.num_particles > 0


def test_native_functors_agree(config):
    states = [ap.AutoPas(config, functors=functor())() for functor in NATIVE_FUNCTORS]
    assert_same_state(*states, atol=1e-7)
    assert states[0].potential_energy == pytest.approx(states[1].potential_energy)


def test_python_functor_equals_native_functor(config):
    functor = PythonLennardJones()
    python = ap.AutoPas(config, functors=functor)()
    native = ap.AutoPas(config, functors=ap.LennardJones())()
    assert functor.calls > 0
    assert_same_state(python, native, atol=1e-8)
    assert python.potential_energy == pytest.approx(native.potential_energy)
    np.testing.assert_allclose(python.history["potential_energy"], native.history["potential_energy"])


@pytest.mark.parametrize("newton3", ["enabled", "disabled"])
def test_energy_is_conserved(newton3):
    """Velocity Verlet has no drift for this step size. This checks integrator, boundaries, and functor in one go."""
    velocities = np.random.default_rng(3).normal(0.0, 1.0, (216, 3))
    crystal = ap.generators.CubeGrid([6, 6, 6], particle_spacing=1.2, bottom_left_corner=[0.1, 0.1, 0.1]).generate()
    config = small_config(iterations=1000, newton3=[newton3], **{"Objects": {}, "data-layout": ["AoS", "SoA"]})
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    sim.add_particles(crystal["positions"], velocities=velocities)
    sim()
    history = sim.history
    total = history["potential_energy"] + history["kinetic_energy"]
    # The first iteration is off by design, see the note on Simulation::run() in md-flexible.
    assert np.ptp(total[10:]) < 1e-2 * abs(total[10:].mean())
    assert history["kinetic_energy"][-1] > 0.0


def test_reflective_boundaries_keep_particles():
    # The wall repels within 2^(1/6) sigma / 2 = 0.56, so the particles have to start further away from it.
    crystal = ap.generators.CubeGrid([6, 6, 6], particle_spacing=1.2, bottom_left_corner=[0.6, 0.6, 0.6])
    config = small_config(
        iterations=200,
        **{"boundary-type": ["reflective"] * 3, "globalForce": [0.0, 0.0, -5.0], "Objects": [crystal]},
    )
    sim = ap.AutoPas(config, functors=ap.LennardJonesHWY())()
    assert sim.num_particles == 216
    assert (sim.positions >= 0.0).all() and (sim.positions <= 7.2).all()


def test_tuning_selects_a_configuration_and_can_be_resumed():
    config = small_config(
        iterations=60,
        **{
            "container": ["LinkedCells", "VerletLists"],
            "traversal": ["lc_c08", "vl_list_iteration"],
            # vl_list_iteration requires this. Without it, the search space would be trivial.
            "newton3": ["enabled", "disabled"],
            "tuning-interval": 30,
            "tuning-samples": 2,
        },
    )
    sim = ap.AutoPas(config, functors=ap.LennardJonesHWY())()
    assert sim.history["tuning"].any()
    assert "pairwise" in sim.current_configuration
    sim(10)
    assert sim.iteration == 70


def test_two_functors_of_a_type_add_up(config):
    class Constant(ap.PairwiseFunctor):
        def aos(self, i, j, newton3):
            i.add_f([1.0, 0.0, 0.0])
            if newton3:
                j.sub_f([1.0, 0.0, 0.0])

    # Only one functor is measured by the tuner. No matter the order of the list, the native one has to run first.
    sim = ap.AutoPas(config, functors=[Constant(), ap.LennardJones()])
    assert [f.name for f in sim.functors] == ["LennardJones", "Constant"]
    sim(2)
    assert np.isnan(sim.potential_energy)
    assert sim.num_particles == 216


def test_triwise_python_functor(config):
    class Counting(ap.TriwiseFunctor):
        triplets = 0

        def aos(self, i, j, k, newton3):
            type(self).triplets += 1

    sim = ap.AutoPas(
        small_config(
            iterations=1,
            **{
                "Objects": [ap.generators.CubeGrid([3, 3, 3], particle_spacing=1.2, bottom_left_corner=[2, 2, 2])],
                "traversal-3b": ["lc_c01"],
                "data-layout-3b": ["AoS"],
                "newton3-3b": ["disabled"],
                "boundary-type": ["none"] * 3,
            },
        ),
        functors=Counting(),
    )
    sim(1)
    assert Counting.triplets > 0


def test_multiple_instances_live_side_by_side(config):
    a = ap.AutoPas(config, functors=ap.LennardJones())
    b = ap.AutoPas(config, functors=ap.LennardJones())
    a(5)
    b(5)
    del b
    a(5)
    assert a.iteration == 10


def test_callback_sees_consistent_state_and_can_stop(config):
    seen = []

    def callback(sim):
        seen.append((sim.iteration, sim.num_particles))
        return sim.iteration < 6

    sim = ap.AutoPas(config, functors=ap.LennardJones())
    sim(100, callback=callback, callback_interval=2)
    assert seen == [(2, 216), (4, 216), (6, 216)]
    assert sim.iteration == 6
    sim()  # Continues to the end of the configuration.
    assert sim.iteration == 20


def test_raising_callback_allows_resuming(config):
    def callback(sim):
        raise KeyError("stop")

    sim = ap.AutoPas(config, functors=ap.LennardJones())
    with pytest.raises(KeyError):
        sim(10, callback=callback)
    assert sim.iteration == 1
    sim(9)
    assert_same_state(sim, ap.AutoPas(config, functors=ap.LennardJones())(10))


def test_nested_run_is_rejected(config):
    sim = ap.AutoPas(config, functors=ap.LennardJones())
    with pytest.raises(RuntimeError, match="already running"):
        sim(1, callback=lambda s: s(1))


def test_python_functor_exception_propagates_and_invalidates_the_state(config):
    class Broken(ap.PairwiseFunctor):
        def aos(self, i, j, newton3):
            raise ZeroDivisionError("from the functor")

    sim = ap.AutoPas(config, functors=Broken())
    with pytest.raises(ZeroDivisionError, match="from the functor"):
        sim(1)
    with pytest.raises(RuntimeError, match="inconsistent"):
        sim(1)
    with pytest.raises(RuntimeError, match="inconsistent"):
        sim.positions


def test_python_functor_without_aos_raises(config):
    class Empty(ap.PairwiseFunctor):
        pass

    with pytest.raises(NotImplementedError, match="override"):
        ap.AutoPas(config, functors=Empty())(1)


def test_cutoff_filter_of_python_functor(config):
    class Distances(ap.PairwiseFunctor):
        largest = 0.0

        def aos(self, i, j, newton3):
            dr = i.r - j.r
            type(self).largest = max(type(self).largest, float(np.sqrt(dr.dot(dr))))

    ap.AutoPas(config, functors=Distances(cutoff=1.5))(1)
    assert 0.0 < Distances.largest <= 1.5
    with pytest.raises(ValueError, match="exceeds"):
        ap.AutoPas(config, functors=Distances(cutoff=CUTOFF + 1))


def test_invalid_input_raises_instead_of_killing_the_interpreter(config, tmp_path):
    with pytest.raises(ValueError, match="Unrecognized option"):
        ap.AutoPas({**config, "no-such-option": 1}, functors=ap.LennardJones())
    with pytest.raises(FileNotFoundError):
        ap.AutoPas(tmp_path / "missing.yaml")
    with pytest.raises(TypeError, match="can not be represented"):
        ap.AutoPas({**config, "cutoff": object()}, functors=ap.LennardJones())
    with pytest.raises(ValueError, match="Python functors need AoS"):
        ap.AutoPas({**config, "data-layout": ["SoA"]}, functors=PythonLennardJones())
    with pytest.raises(ValueError, match="one native functor"):
        ap.AutoPas(config, functors=[ap.LennardJones(), ap.LennardJonesHWY()])
    with pytest.raises(ValueError, match="unknown type"):
        sites = {0: {"mass": 1.0, "epsilon": 1.0, "sigma": 1.0}}
        objects = [ap.generators.CubeGrid([2, 2, 2], particle_spacing=1.2, bottom_left_corner=[1, 1, 1], type_id=3)]
        ap.AutoPas({**config, "Sites": sites, "Objects": objects}, functors=ap.LennardJones())


def test_functor_can_not_be_shared(config):
    functor = ap.LennardJones()
    ap.AutoPas(config, functors=functor)
    with pytest.raises(ValueError, match="already used"):
        ap.AutoPas(config, functors=functor)


def test_add_particles(config):
    sim = ap.AutoPas({**config, "Objects": {}}, functors=ap.LennardJones())
    assert sim.num_particles == 0
    sim.add_particles([[1.0, 1.0, 1.0], [2.2, 1.0, 1.0]], velocities=[0.5, 0.0, 0.0])
    sim.add_particles(np.array([[4.0, 4.0, 4.0]]), ids=[7])
    assert sim.num_particles == 3
    np.testing.assert_array_equal(sim.ids, [0, 1, 7])
    np.testing.assert_allclose(sim.velocities[1], [0.5, 0.0, 0.0])
    with pytest.raises(ValueError, match="outside"):
        sim.add_particles([[1.0, 1.0, 1.0], [100.0, 1.0, 1.0]])
    assert sim.num_particles == 3
    with pytest.raises(ValueError, match="unknown type"):
        sim.add_particles([[1.0, 1.0, 1.0]], type_ids=1)
    sim(5)
    assert sim.iteration == 5


def test_history(config):
    sim = ap.AutoPas(config, functors=ap.LennardJones(), history_interval=5)()
    history = sim.history
    np.testing.assert_array_equal(history["iteration"], [0, 5, 10, 15])
    assert history["tuning"].dtype == bool
    assert sim.temperature == pytest.approx(sim.kinetic_energy / (1.5 * sim.num_particles))
    assert len(ap.AutoPas(config, functors=ap.LennardJones(), history_interval=0)().history["iteration"]) == 0


def periodic_reference(positions, box_length, cutoff=CUTOFF):
    """Truncated and shifted Lennard-Jones for sigma = epsilon = 1 by brute force with the minimum image convention."""
    dr = positions[:, None] - positions[None]
    dr -= box_length * np.round(dr / box_length)
    r2 = (dr * dr).sum(-1)
    np.fill_diagonal(r2, np.inf)
    within_cutoff = r2 < cutoff**2
    inv_r6 = 1.0 / r2**3
    energy = 0.5 * ((4.0 * inv_r6 * (inv_r6 - 1.0) - 4.0 * (cutoff**-12 - cutoff**-6)) * within_cutoff).sum()
    forces = ((24.0 * (2.0 * inv_r6**2 - inv_r6) / r2 * within_cutoff)[..., None] * dr).sum(1)
    return energy, forces


CONFIGURATIONS = [
    ("LinkedCells", "lc_c08", "enabled"),
    ("LinkedCells", "lc_c01", "disabled"),
    ("VerletLists", "vl_list_iteration", "disabled"),
    ("VerletListsCells", "vlc_c18", "enabled"),
    ("VerletClusterLists", "vcl_c06", "enabled"),
    ("DirectSum", "ds_sequential", "disabled"),
]


def evaluate_once(functor, container, traversal, newton3, data_layout):
    """One force evaluation of a perturbed crystal. The time step is non-zero because halos are skipped otherwise."""
    box_length = 7.2
    positions = np.mod(
        np.mgrid[0:6, 0:6, 0:6].reshape(3, -1).T * 1.2 + 0.1 + np.random.default_rng(1).normal(0.0, 0.08, (216, 3)),
        box_length,
    )
    config = small_config(
        iterations=1,
        deltaT=1e-12,
        container=[container],
        traversal=[traversal],
        newton3=[newton3],
        Objects={},
        **{"data-layout": [data_layout]},
    )
    sim = ap.AutoPas(config, functors=functor)
    sim.add_particles(positions)
    sim(1)
    return sim, *periodic_reference(sim.positions, box_length)


@pytest.mark.parametrize("data_layout", ["AoS", "SoA"])
@pytest.mark.parametrize("container, traversal, newton3", CONFIGURATIONS)
def test_lennard_jones_against_brute_force(container, traversal, newton3, data_layout):
    sim, energy, forces = evaluate_once(ap.LennardJones(), container, traversal, newton3, data_layout)
    np.testing.assert_allclose(sim.forces, forces, atol=1e-9)
    assert sim.potential_energy == pytest.approx(energy)


@pytest.mark.filterwarnings("ignore:LennardJonesHWY reports")
@pytest.mark.parametrize("data_layout", ["AoS", "SoA"])
@pytest.mark.parametrize("container, traversal, newton3", CONFIGURATIONS)
def test_highway_forces_against_brute_force(container, traversal, newton3, data_layout):
    sim, _, forces = evaluate_once(ap.LennardJonesHWY(), container, traversal, newton3, data_layout)
    np.testing.assert_allclose(sim.forces, forces, atol=1e-9)


@pytest.mark.filterwarnings("ignore:LennardJonesHWY reports")
@pytest.mark.parametrize(
    "data_layout",
    [
        "AoS",
        pytest.param(
            "SoA",
            marks=pytest.mark.xfail(
                strict=True,
                reason="LJFunctorHWY: The SoA kernels mistake non-dummy for owned, so halo pairs count fully.",
            ),
        ),
    ],
)
def test_highway_energy_against_brute_force(data_layout):
    sim, energy, _ = evaluate_once(ap.LennardJonesHWY(), "LinkedCells", "lc_c08", "enabled", data_layout)
    assert sim.potential_energy == pytest.approx(energy)


def test_warning_about_highway_energies(config):
    with pytest.warns(UserWarning, match="LennardJonesHWY reports"):
        ap.AutoPas({**config, "data-layout": ["AoS", "SoA"]}, functors=ap.LennardJonesHWY())


def test_python_functors_run_on_one_thread_and_leave_others_untouched(config):
    threads = set()

    class Threads(ap.PairwiseFunctor):
        def aos(self, i, j, newton3):
            threads.add(threading.get_ident())

    before = ap.max_threads()
    ap.AutoPas(config, functors=Threads())(2)
    assert len(threads) == 1
    assert ap.max_threads() == before
