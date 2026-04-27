#! /usr/bin/env python3
from pathlib import Path
from string import Template


SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE_FILE = SCRIPT_DIR / "rayleighTaylor.yaml"
OUTPUT_ROOT = SCRIPT_DIR / "generated_inputs"

TRAVERSALS = ["hgrid_matching", "hgrid_block4", "hgrid_block8", "lc_c08"]
CELL_SIZES = [0.5, 1.0]

TRAVERSAL_TO_CONTAINER = {
    "hgrid_matching": "HierarchicalGridMatching",
    "hgrid_block4": "HierarchicalGrid",
    "hgrid_block8": "HierarchicalGrid",
    "lc_c08": "LinkedCells",
}


def format_value(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


input_template = Template(TEMPLATE_FILE.read_text(encoding="utf-8"))
OUTPUT_ROOT.mkdir(exist_ok=True)

# Minimal folder layout:
# generated_inputs/traversal_<...>/cellSize_<...>/input.yaml
for traversal in TRAVERSALS:
    container = TRAVERSAL_TO_CONTAINER[traversal]
    traversal_dir = OUTPUT_ROOT / f"traversal_{traversal}"
    traversal_dir.mkdir(exist_ok=True)

    for cell_size in CELL_SIZES:
        cell_size_dir = traversal_dir / f"cellSize_{format_value(cell_size)}"
        cell_size_dir.mkdir(exist_ok=True)

        substitutions = {
            "container": container,
            "traversal": traversal,
            "cellSize": f"{cell_size}",
        }

        (cell_size_dir / "input.yaml").write_text(input_template.substitute(substitutions), encoding="utf-8")