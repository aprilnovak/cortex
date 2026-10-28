# TEA — Techno-Economic Analysis cost model

A manufacturing cost model for fusion reactor components, based on the
Swift Design-for-Manufacture (DFM) relative-cost-coefficient
method: a component's fabrication cost per process is `Rc * Pc`, where `Pc`
is a base processing-cost curve and `Rc = Cmp * Cc * Cs * max(Ct, Cf)`
combines material/process compatibility (`Cmp`) with geometry-driven
complexity/size/tolerance/finish coefficients (`Cc`/`Cs`/`Ct`/`Cf`).

There are two ways to run it:

* **`tea_api`** — a programmatic Python API. Import it in your own script or
  notebook to evaluate one material or screen a batch of candidates.
* **The web GUI** — a FastAPI backend + single-page frontend for interactive,
  point-and-click use.

Both sit on the same underlying modules (`tea1.py`, `material_def.py`,
`process_def.py`, `component_def.py`, `cost_variables.py`), so results from
one match the other exactly.

## Install

`tea/` is a pip-installable package. From this directory:

```bash
pip install -e .
```

The core cost model has no third-party dependencies, so that's all you need
for `tea_api`. `-e` ("editable") installs it in place, pointing back at this
checkout rather than copying files.

Optional extras, if you need them:

```bash
pip install -e ".[webapp]"    # fastapi, uvicorn, matplotlib, pydantic — to run the web GUI
pip install -e ".[notebook]"  # pandas, ipykernel — for the example notebook's optional cells
```

## Using `tea_api`

Start with [`examples/tea_api_example.ipynb`](examples/tea_api_example.ipynb) —
a guided walkthrough covering single-material evaluation, custom material
compositions, blending, custom processes, geometry selection, batch
screening (a list or a whole directory of candidate `.py` files), and a
cheat sheet of every public function.

Minimal example:

```python
import tea_api

result = tea_api.evaluate(
    "EROFER97",              # material name, looked up in materials_database/
    ["CNC", "Hot Rolling"],  # process names, looked up in processes_database/
    volume_mm3=3_000_000,
    production_qty=100,
)
print(result["summary"]["Total cost"])
```

## Running the web GUI

```bash
pip install -e ".[webapp]"
cd webapp
uvicorn server:app --reload
```

Then open `http://127.0.0.1:8000` — it's a single self-contained
`static/index.html` (Vue via CDN, no build step) talking to the FastAPI
backend's `/api/*` routes, which wrap the same `tea_api` functions used
above.

## Directory layout

* `materials_database/`, `processes_database/`, `geometries_database/` —
  the material/process/geometry `.py` definition files `tea_api` and the
  web GUI look names up in by default. Each is a plain module with
  module-level variables (`name`, `density`, `composition`, `Cmp_map`, ...
  — see any existing file for the format). Point `tea_api` at a different
  directory instead with the `material_dir=`/`process_dir=`/`geometry_dir=`
  keyword arguments.
* `inputs_materials/`, `inputs_processes/`, `inputs_components/` — example/
  template input files in the same format, used by the archived interactive
  workflow below.
* `inputs_materials_machiningdoctor/` — a large set of reference material
  property files sourced from the MachiningDoctor database.
* `archived_modules/` — the original interactive, `input()`-prompt-driven
  scripts (`tea_script.py` and friends) this project has moved on from, kept
  for reference.
* `examples/` — the `tea_api` walkthrough notebook.
* `webapp/` — the FastAPI backend (`server.py`) and static frontend
  (`static/index.html`) for the GUI.
