"""
FastAPI backend for the TEA web GUI.

Wraps the existing CLI-oriented modules (material_def, process_def, tea1) in
a small JSON API and serves the Vue frontend from ./static. Run from this
directory with:

    uvicorn server:app --reload
"""

import base64
import csv
import io
import os
import sys
from pathlib import Path
from typing import Dict, List, Literal, Optional

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel, Field

# ----------------------------------------------------------------------
# Make the cortex_tea package importable without a pip install, and run from
# tea/ so relative user-input directories like "inputs_materials" resolve.

TEA_DIR = Path(__file__).resolve().parent.parent
os.chdir(TEA_DIR)
sys.path.insert(0, str(TEA_DIR))

from cortex_tea import material_def
from cortex_tea import process_def
from cortex_tea import component_def
from cortex_tea import cost_variables
from cortex_tea import tea_api
from cortex_tea import plotting
from cortex_tea.process_def import (
    COST_LEVELS,
    TOOLING_COST_BY_LEVEL,
    EQUIPMENT_COST_BY_LEVEL,
    PROCESSING_TIME_BY_LEVEL,
    EQUIPMENT_UTILIZATION_SEC,
    EQUIPMENT_ESCALATION,
)
from cortex_tea.tea1 import Component

STATIC_DIR = Path(__file__).resolve().parent / "static"

app = FastAPI(title="TEA Web GUI")

# ----------------------------------------------------------------------
# Helpers
#
# Material/process/geometry lookup, resolution, and serialization now
# live in tea_api.py (shared with the programmatic researcher-facing API);
# this module just translates tea_api.TeaError/TeaNotFoundError into
# HTTPException.


# ----------------------------------------------------------------------
# Categorical cost levels for custom processes now live in process_def.py
# (COST_LEVELS, *_BY_LEVEL tables, compute_alphaT_beta) so the same logic
# is shared with the process database files and any future CLI use.


# ----------------------------------------------------------------------
# Request/response models


class CompositionElement(BaseModel):
    wt: float = Field(gt=0)
    type: Literal["alloy", "residual"]


class MaterialInput(BaseModel):
    mode: Literal["existing", "custom", "blend"]
    name: str
    material_dir: Optional[str] = None
    density: Optional[float] = None
    composition: Optional[Dict[str, CompositionElement]] = None
    remainder_element: Optional[str] = None
    material_a_name: Optional[str] = None
    material_b_name: Optional[str] = None
    blend_fraction_a: Optional[float] = Field(default=None, gt=0, lt=1)
    blend_basis: Literal["volume", "mass"] = "volume"


CostLevel = Literal["Low", "Low-Medium", "Medium", "Medium-High", "High", "High-Very High", "Very High"]


class ProcessInput(BaseModel):
    mode: Literal["existing", "custom"]
    name: str
    process_dir: Optional[str] = None
    tooling_level: Optional[CostLevel] = None
    equipment_level: Optional[CostLevel] = None
    time_level: Optional[CostLevel] = None
    description: Optional[str] = ""
    Cc: float = Field(default=1.0, gt=0)
    Cs: float = Field(default=1.0, gt=0)
    Ct: float = Field(default=1.0, gt=0)
    Cf: float = Field(default=1.0, gt=0)
    Wc: float = Field(default=1.0, gt=0)
    Cmp: Optional[float] = Field(default=None, gt=0)


class GeometryInput(BaseModel):
    name: str
    geometry_dir: Optional[str] = None


class CalculateRequest(BaseModel):
    component_name: Optional[str] = None
    volume_mm3: float = Field(gt=0)
    production_qty: int = Field(gt=0)
    material: MaterialInput
    processes: List[ProcessInput]
    geometry: Optional[GeometryInput] = None


class ScreenDirectoryRequest(BaseModel):
    component_name: Optional[str] = None
    volume_mm3: float = Field(gt=0)
    production_qty: int = Field(gt=0)
    material_dir: str
    processes: List[ProcessInput]
    geometry: Optional[GeometryInput] = None


# ----------------------------------------------------------------------
# Routes


@app.get("/api/materials")
def get_materials(material_dir: Optional[str] = None):
    try:
        lookup = tea_api._load_material_lookup(material_dir)
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))
    return [tea_api.serialize_material(m) for m in sorted(lookup.values(), key=lambda m: m.name)]


@app.get("/api/processes")
def get_processes(process_dir: Optional[str] = None):
    try:
        lookup = tea_api._load_process_lookup(process_dir)
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))
    return [tea_api.serialize_process(p) for p in sorted(lookup.values(), key=lambda p: p.name)]


@app.get("/api/element-costs")
def get_element_costs():
    return {el: mc.cost for el, mc in cost_variables.material_cost_database.items()}


@app.get("/api/geometries")
def get_geometries(geometry_dir: Optional[str] = None):
    try:
        lookup = tea_api._load_geometry_lookup(geometry_dir)
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))
    return [tea_api.serialize_geometry(g) for g in sorted(lookup.values(), key=lambda g: g.name)]


@app.get("/api/cost-levels")
def get_cost_levels():
    return {
        "levels": COST_LEVELS,
        "tooling_cost": TOOLING_COST_BY_LEVEL,
        "equipment_cost": EQUIPMENT_COST_BY_LEVEL,
        "processing_time": PROCESSING_TIME_BY_LEVEL,
        "equipment_utilization_sec": EQUIPMENT_UTILIZATION_SEC,
        "equipment_escalation": EQUIPMENT_ESCALATION,
    }


def _material_spec(material_input: "MaterialInput") -> "tea_api.MaterialSpec":
    """Translate a MaterialInput (any of the three modes) into the plain
    str/dict spec tea_api.resolve_material() expects. Shared by
    /api/calculate and /api/material-preview so both resolve a draft
    material identically."""
    if material_input.mode == "existing":
        return material_input.name

    if material_input.mode == "custom":
        composition = {el: {"wt": c.wt, "type": c.type} for el, c in (material_input.composition or {}).items()}
        return {
            "name": material_input.name,
            "density": material_input.density,
            "composition": composition,
            "remainder_element": material_input.remainder_element,
        }

    # blend
    return {
        "name": material_input.name,
        "material_a": material_input.material_a_name,
        "material_b": material_input.material_b_name,
        "fraction_a": material_input.blend_fraction_a,
        "basis": material_input.blend_basis,
    }


def _process_spec(process_input: "ProcessInput") -> "tea_api.ProcessSpec":
    """Translate a ProcessInput into the plain dict spec
    tea_api.resolve_process() expects, carrying along its per-process
    Cc/Cs/Ct/Cf/Wc/Cmp overrides for build_component(). For "existing" mode
    the categorical-level keys are omitted entirely (not just set to None)
    so resolve_process() takes the name-lookup branch rather than the
    custom-build branch."""
    overrides = {
        "Cc": process_input.Cc, "Cs": process_input.Cs, "Ct": process_input.Ct,
        "Cf": process_input.Cf, "Wc": process_input.Wc, "Cmp": process_input.Cmp,
    }
    if process_input.mode == "existing":
        return {"name": process_input.name, **overrides}

    return {
        "name": process_input.name,
        "tooling_level": process_input.tooling_level,
        "equipment_level": process_input.equipment_level,
        "time_level": process_input.time_level,
        "description": process_input.description or "",
        **overrides,
    }


@app.post("/api/material-preview")
def material_preview(material_input: MaterialInput):
    """Authoritative preview of a draft material's density and composition
    (existing/custom/blend) - the frontend only displays this, it does not
    recompute any of it. No chart here: charts render once, from the final
    resolved material, as part of /api/calculate's response (Results)."""
    try:
        material = tea_api.resolve_material(_material_spec(material_input), material_dir=material_input.material_dir)
    except tea_api.TeaNotFoundError as e:
        raise HTTPException(404, str(e))
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))

    composition, warnings = tea_api.material_cost_breakdown(material)
    return {
        "name": material.name,
        "density": material.density,
        "cost": material.cost,
        "composition": composition,
        "warnings": warnings,
    }


@app.post("/api/process-preview")
def process_preview(process_input: ProcessInput):
    """Authoritative preview of a draft custom process's alphaT/beta - the
    frontend only displays this, it does not recompute compute_alphaT_beta()
    itself."""
    try:
        proc = tea_api.resolve_process(_process_spec(process_input), process_dir=process_input.process_dir)
    except tea_api.TeaNotFoundError as e:
        raise HTTPException(404, str(e))
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))
    return {"name": proc.name, "alphaT": proc.alphaT, "beta": proc.beta}


@app.post("/api/calculate")
def calculate(req: CalculateRequest):
    if not req.processes:
        raise HTTPException(400, "At least one process must be specified.")

    try:
        material = tea_api.resolve_material(_material_spec(req.material), material_dir=req.material.material_dir)

        geometry = None
        if req.geometry is not None:
            geometry = tea_api.resolve_geometry(req.geometry.name, geometry_dir=req.geometry.geometry_dir)

        process_specs = [_process_spec(p) for p in req.processes]
        processes = [
            tea_api.resolve_process(spec, process_dir=p.process_dir)
            for spec, p in zip(process_specs, req.processes)
        ]

        component = tea_api.build_component(
            req.component_name or "Component", req.volume_mm3, req.production_qty,
            material, processes, process_specs=process_specs, geometry=geometry,
        )
    except tea_api.TeaNotFoundError as e:
        raise HTTPException(404, str(e))
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))

    try:
        summary = component.manufacturing_cost()
    except Exception as exc:
        raise HTTPException(400, str(exc))

    _, warnings = tea_api.material_cost_breakdown(material)

    chart_png = plotting.cost_breakdown_png(component)
    curves_png = plotting.process_cost_curves_png(component)

    return {
        "summary": summary,
        "chart_png_base64": base64.b64encode(chart_png).decode("ascii"),
        "process_cost_curve_chart_png_base64": base64.b64encode(curves_png).decode("ascii"),
        "warnings": warnings,
    }


@app.post("/api/screen-directory")
def screen_directory_route(req: ScreenDirectoryRequest):
    """Run tea_api.screen_directory() over every material file in
    req.material_dir against the same processes/geometry/volume/production
    quantity a normal /api/calculate call would use, and return the results
    as a downloadable CSV (one row per material)."""
    if not req.processes:
        raise HTTPException(400, "At least one process must be specified.")

    process_specs = [_process_spec(p) for p in req.processes]

    try:
        results = tea_api.screen_directory(
            req.material_dir, process_specs, req.volume_mm3, req.production_qty,
            geometry=req.geometry.name if req.geometry else None,
            geometry_dir=req.geometry.geometry_dir if req.geometry else None,
            component_name=req.component_name or "Component",
        )
    except tea_api.TeaNotFoundError as e:
        raise HTTPException(404, str(e))
    except tea_api.TeaError as e:
        raise HTTPException(400, str(e))

    process_names = [p.name for p in req.processes]

    # Union of elements across every successful candidate, in first-seen
    # order, so the elemental cost breakdown lines up in one fixed set of
    # columns even though different candidates can have different
    # compositions.
    element_names = []
    seen_elements = set()
    for row in results:
        if row["ok"]:
            for elem_row in row["result"]["material_composition"]:
                if elem_row["element"] not in seen_elements:
                    seen_elements.add(elem_row["element"])
                    element_names.append(elem_row["element"])

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(
        ["Material", "Status", "Error", "Total Cost ($)", "Unit Cost ($/kg)", "Mass (kg)",
         "Material Cost ($)", "Processing Cost ($)"]
        + [f"{name} Cost ($)" for name in process_names]
        + [f"{el} Cost ($/kg)" for el in element_names]
    )
    for row in results:
        if row["ok"]:
            s = row["result"]["summary"]
            breakdown = s["Cost Breakdown"]
            elem_costs = {r["element"]: r["dollar_per_kg"] for r in row["result"]["material_composition"]}
            writer.writerow(
                [row["material"], "OK", "", f"{s['Total cost']:.2f}", f"{s['Unit cost']:.2f}",
                 f"{s['Mass']:.4f}", f"{s['Material cost']:.2f}", f"{s['Processing cost']:.2f}"]
                + [f"{breakdown[name]:.2f}" if name in breakdown else "" for name in process_names]
                + [f"{elem_costs[el]:.4f}" if el in elem_costs else "" for el in element_names]
            )
        else:
            writer.writerow(
                [row["material"], "FAILED", row["error"], "", "", "", "", ""]
                + [""] * len(process_names)
                + [""] * len(element_names)
            )

    return Response(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=tea_screen_results.csv"},
    )


# ----------------------------------------------------------------------
# Static frontend


@app.get("/")
def index():
    return FileResponse(STATIC_DIR / "index.html")


app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
