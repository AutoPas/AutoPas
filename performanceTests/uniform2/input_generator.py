#! /usr/bin/env python3
from pathlib import Path
from string import Template


# Uses template_input.yaml to create input files that emulate a simulation
# domain distributed uniformly with two particle types.
#
# This sweep varies:
# - container / traversal combinations
# - sigma ratio: sigma0 / sigmaMax
# - particle-count ratio: n1 / n0
# - cell size
#
# For each scenario, NUM_REPEATS repeat runs are generated.

SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_FILE = SCRIPT_DIR / "template_input.yaml"
OUTPUT_ROOT = SCRIPT_DIR / "generated_inputs"

NUM_REPEATS = 3
TOTAL_PARTICLES = 1_000_000
SIGMAMAX = 0.8

SIGMA_RATIOS = [0.5,0.6, 0.7, 0.8, 0.9]
COUNT_RATIOS = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
CELL_SIZES = [0.5, 1.0]

# Container -> traversals to sweep.
CONTAINER_TRAVERSALS = {
    "HierarchicalGridMatching": ["hgrid_matching"],
    "HierarchicalGrid": ["hgrid_block4", "hgrid_block8"],
    "LinkedCells": ["lc_c08"],
}

def format_ratio(value: float) -> str:
    """Convert ratio to a filesystem-friendly value."""
    return f"{value:.2f}".replace(".", "p")


def split_particles(total_particles, count_ratio):
    """Split total particles into n0 and n1 for the given n1 / n0 ratio."""
    n1 = max(1, int(round(total_particles / (1.0 + count_ratio))))
    n0 = total_particles - n1
    return n0, n1

CONTAINER_TRAVERSALS = {
    "HierarchicalGridMatching": ["hgrid_matching"],}

input_template = Template(TEMPLATE_FILE.read_text(encoding="utf-8"))

OUTPUT_ROOT.mkdir(exist_ok=True)

# Create directory structure:
# generated_inputs/container_<container>/traversal_<traversal>/sigmaRatio_<...>/cellSize_<...>/countRatio_<...>/run_<...>
for container, traversals in CONTAINER_TRAVERSALS.items():
    container_dir = OUTPUT_ROOT / f"container_{container}"
    container_dir.mkdir(exist_ok=True)

    for traversal in traversals:
        traversal_dir = container_dir / f"traversal_{traversal}"
        traversal_dir.mkdir(exist_ok=True)

        for sigma_ratio in SIGMA_RATIOS:
            sigma0 = SIGMAMAX * sigma_ratio
            sigma1 = SIGMAMAX
            cutoff0 = 3.0 * sigma0
            cutoff1 = 3.0 * sigma1
            max_cutoff = max(cutoff0, cutoff1)

            sigma_dir = traversal_dir / f"sigmaRatio_{format_ratio(sigma_ratio)}"
            sigma_dir.mkdir(exist_ok=True)

            for cell_size in CELL_SIZES:
                cell_size_dir = sigma_dir / f"cellSize_{format_ratio(cell_size)}"
                cell_size_dir.mkdir(exist_ok=True)

                for count_ratio in COUNT_RATIOS:
                    num_particles0, num_particles1 = split_particles(TOTAL_PARTICLES, count_ratio)
                    count_dir = cell_size_dir / f"countRatio_{format_ratio(count_ratio)}"
                    count_dir.mkdir(exist_ok=True)

                    for run in range(NUM_REPEATS):
                        run_dir = count_dir / f"run_{run}"
                        run_dir.mkdir(exist_ok=True)

                        substitutions = {
                            "sigma0": f"{sigma0:.6f}",
                            "sigma1": f"{sigma1:.6f}",
                            "cutoff0": f"{cutoff0:.6f}",
                            "cutoff1": f"{cutoff1:.6f}",
                            "maxCutoff": f"{max_cutoff:.6f}",
                            "numParticles0": str(num_particles0),
                            "numParticles1": str(num_particles1),
                            "container": container,
                            "traversal": traversal,
                            "cellSize": f"{cell_size}",
                        }

                        (run_dir / "input.yaml").write_text(
                            input_template.substitute(substitutions), encoding="utf-8"
                        )

