#! /usr/bin/env python3
from pathlib import Path
from string import Template
import numpy as np


# Uses template_input.yaml to create input files that emulate a simulation
# domain distributed uniformly with two particle types.
#
# This sweep varies:
# - sigma ratio: sigma1 / sigma0
# - data layout: AoS or SoA
#
# The particle-count ratio n1 / n0 is fixed to 2.0.
# For each scenario, 3 repeat runs are generated.

SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_FILE = SCRIPT_DIR / "template_input.yaml"
OUTPUT_ROOT = SCRIPT_DIR / "generated_inputs_rangeCheck"

NUM_REPEATS = 3
TOTAL_PARTICLES = 120_000
SIGMAMAX = 1.0

SIGMA_RATIOS = np.linspace(0.05, 0.55, 6)
COUNT_RATIOS = [2.0]
DATA_LAYOUTS = ["AoS", "SoA"]


def format_ratio(value: float) -> str:
    """Convert ratio to a filesystem-friendly value."""
    return f"{value:.2f}".replace(".", "p")


def split_particles(total_particles, count_ratio):
    """Split total particles into n0 and n1 for the given n1 / n0 ratio."""
    n1 = max(1, int(round(total_particles / (1.0 + count_ratio))))
    n0 = total_particles - n1
    return n0, n1


input_template = Template(TEMPLATE_FILE.read_text(encoding="utf-8"))

OUTPUT_ROOT.mkdir(exist_ok=True)

# Create directory structure:
# generated_inputs_rangeCheck/sigmaRatio_<sigma1_over_sigma0>/countRatio_<n1_over_n0>/dataLayout_<AoS|SoA>/run_<run>
for sigma_ratio in SIGMA_RATIOS:
    sigma0 = SIGMAMAX * sigma_ratio
    sigma1 = SIGMAMAX
    cutoff0 = 3.0 * sigma0
    cutoff1 = 3.0 * sigma1
    max_cutoff = max(cutoff0, cutoff1)

    sigma_dir = OUTPUT_ROOT / f"sigmaRatio_{format_ratio(sigma_ratio)}"
    sigma_dir.mkdir(exist_ok=True)

    for count_ratio in COUNT_RATIOS:
        num_particles0, num_particles1 = split_particles(TOTAL_PARTICLES, count_ratio)
        count_dir = sigma_dir / f"countRatio_{format_ratio(count_ratio)}"
        count_dir.mkdir(exist_ok=True)

        for data_layout in DATA_LAYOUTS:
            layout_dir = count_dir / f"dataLayout_{data_layout}"
            layout_dir.mkdir(exist_ok=True)

            for run in range(NUM_REPEATS):
                run_dir = layout_dir / f"run_{run}"
                run_dir.mkdir(exist_ok=True)

                substitutions = {
                    "sigma0": f"{sigma0:.6f}",
                    "sigma1": f"{sigma1:.6f}",
                    "cutoff0": f"{cutoff0:.6f}",
                    "cutoff1": f"{cutoff1:.6f}",
                    "maxCutoff": f"{max_cutoff:.6f}",
                    "numParticles0": str(num_particles0),
                    "numParticles1": str(num_particles1),
                    "dataLayout": data_layout,
                }

                (run_dir / "input.yaml").write_text(
                    input_template.substitute(substitutions), encoding="utf-8"
                )
