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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import ConnectionPatch

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel, Field

# ----------------------------------------------------------------------
# Wire up imports to the sibling tea/*.py modules (they use bare, sibling
# imports like `from process_def import Process`, so tea/ must be on
# sys.path and also be the cwd for the relative database/input directories
# like "materials_database" to resolve the way tea_script.py expects).

TEA_DIR = Path(__file__).resolve().parent.parent
os.chdir(TEA_DIR)
sys.path.insert(0, str(TEA_DIR))

import material_def
import process_def
import component_def
import cost_variables
import tea_api
from process_def import (
    COST_LEVELS,
    TOOLING_COST_BY_LEVEL,
    EQUIPMENT_COST_BY_LEVEL,
    PROCESSING_TIME_BY_LEVEL,
    EQUIPMENT_UTILIZATION_SEC,
    EQUIPMENT_ESCALATION,
)
from tea1 import Component

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
# Chart rendering
#
# Styling per the dataviz skill: fixed-order categorical palette (never
# cycled), hairline recessive gridlines, muted axis ink, a legend whenever
# 2+ series/segments are on screen, direct value labels where they fit.

_INK_PRIMARY = "#0b0b0b"
_INK_SECONDARY = "#52514e"
_INK_MUTED = "#898781"
_GRIDLINE = "#e1e0d9"
_BASELINE = "#c3c2b7"
_SURFACE = "#fcfcfb"
_BAR_COLOR = "#2a78d6"

# Fixed-order categorical palette (light mode), per the dataviz skill.
# Slot 1 (blue) is reused as _BAR_COLOR for continuity with the prior chart.
_CATEGORICAL_PALETTE = [
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
]
_OTHER_COLOR = _INK_MUTED  # fold-in color for segments beyond the palette

# Sequential single-hue ramp (blue), light -> dark, per the dataviz skill.
# Slot 1 of _CATEGORICAL_PALETTE is step 450 - reused here as the first
# (largest-share) step so the zoom panel reads as "inside" that segment.
_BLUE_RAMP = ["#2a78d6", "#5598e7", "#6da7ec", "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"]


def _label_ink(hex_color: str) -> str:
    """White or ink text, whichever contrasts with the given fill."""
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "#ffffff" if luminance < 140 else _INK_PRIMARY


def _draw_stacked_segments(ax, labels: list, values: list, colors: list, bar_width: float, value_fmt) -> tuple:
    """Draw one stacked bar of (label, value, color) onto ax. Returns the
    segment text artists (for later fit-checking) and each segment's
    baseline y (so a caller can find e.g. a specific segment's top/bottom).
    No gap between segments - the stack's total height is exactly
    sum(values), so a 100%-composition bar tops out at exactly 100%."""
    bottom = 0.0
    texts, bottoms = [], []
    for label, value, color in zip(labels, values, colors):
        ax.bar(0, value, width=bar_width, bottom=bottom, color=color, edgecolor="none", linewidth=0, zorder=3, label=label)
        text = ax.text(
            0, bottom + value / 2, value_fmt(value),
            ha="center", va="center", fontsize=9, color=_label_ink(color), zorder=4,
        )
        texts.append(text)
        bottoms.append(bottom)
        bottom += value
    return texts, bottoms


def _drop_oversized_labels(fig, ax, bar_width: float, texts: list, values: list) -> None:
    """A label that won't fit doesn't get clipped - measure it against the
    segment's rendered width AND height and drop it if it overflows either
    way (the legend still carries the value for segments too small to hold
    a label)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bar_left_px = ax.transData.transform((-bar_width / 2, 0))[0]
    bar_right_px = ax.transData.transform((bar_width / 2, 0))[0]
    available_width_px = (bar_right_px - bar_left_px) - 8  # padding on both sides
    for text, value in zip(texts, values):
        segment_height_px = ax.transData.transform((0, value))[1] - ax.transData.transform((0, 0))[1]
        bbox = text.get_window_extent(renderer=renderer)
        if bbox.width > available_width_px or bbox.height > segment_height_px - 4:
            text.remove()


def _style_stacked_axis(ax, title: str, ylabel: str, value_fmt) -> None:
    if title:
        ax.set_title(title, color=_INK_PRIMARY, fontsize=12, pad=12)
    ax.set_ylabel(ylabel, color=_INK_SECONDARY, fontsize=10)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks([])
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: value_fmt(v)))
    ax.tick_params(axis="y", colors=_INK_MUTED, labelsize=9)
    ax.yaxis.grid(True, color=_GRIDLINE, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    for spine_name in ("top", "right", "left", "bottom"):
        ax.spines[spine_name].set_visible(False)


def render_cost_chart(summary: dict, material_fraction_rows: list) -> str:
    """
    Stacked bar of Material + each process's cost, with a "zoom" panel next
    to it showing the Material segment's cost broken down by element -
    shaded in the same blue as the Material segment (light -> dark, a
    sequential one-hue ramp) and linked to it with connector lines, so it
    reads as a detail view of that one segment (matplotlib's "bar of pie"
    pattern, bar-to-bar instead of pie-to-bar).
    """
    labels = ["Material"] + [p["Process"] for p in summary["Processes"]]
    values = [summary["Cost Breakdown"][label] for label in labels]
    colors = [_CATEGORICAL_PALETTE[i] if i < len(_CATEGORICAL_PALETTE) else _OTHER_COLOR for i in range(len(labels))]

    visible = [r for r in material_fraction_rows if r["fraction"] > 0]
    element_labels = [r["element"] for r in visible]
    element_values = [r["fraction"] * 100 for r in visible]
    element_colors = [_BLUE_RAMP[i] if i < len(_BLUE_RAMP) else _OTHER_COLOR for i in range(len(element_labels))]

    bar_width = 1.2
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 5.4), dpi=150, gridspec_kw={"width_ratios": [1.3, 1]})
    fig.patch.set_facecolor(_SURFACE)
    ax1.set_facecolor(_SURFACE)
    ax2.set_facecolor(_SURFACE)

    texts1, bottoms1 = _draw_stacked_segments(ax1, labels, values, colors, bar_width, value_fmt=lambda v: f"${v:,.2f}")
    material_bottom, material_top = bottoms1[0], bottoms1[0] + values[0]

    texts2, bottoms2 = _draw_stacked_segments(
        ax2, element_labels, element_values, element_colors, bar_width, value_fmt=lambda v: f"{v:.1f}%",
    )
    zoom_top = (bottoms2[-1] + element_values[-1]) if element_values else 0.0

    # No titles here - the page supplies one caption below both panels.
    _style_stacked_axis(ax1, None, "Cost ($)", lambda v: f"${v:,.0f}")
    _style_stacked_axis(ax2, None, "Share of material cost (%)", lambda v: f"{v:.0f}%")

    # Connector lines linking the Material segment (ax1) to the zoom bar (ax2).
    for y1, y2 in ((material_top, zoom_top), (material_bottom, 0.0)):
        con = ConnectionPatch(
            xyA=(bar_width / 2, y1), coordsA=ax1.transData,
            xyB=(-bar_width / 2, y2), coordsB=ax2.transData,
            color=_BASELINE, linewidth=1, linestyle="--", zorder=1,
        )
        fig.add_artist(con)

    legend1 = ax1.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.07), frameon=False, ncol=2,
        fontsize=8, labelcolor=_INK_SECONDARY, handlelength=1.2, handleheight=1.2,
    )
    legend2 = ax2.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.07), frameon=False, ncol=2,
        fontsize=8, labelcolor=_INK_SECONDARY, handlelength=1.2, handleheight=1.2,
    )

    fig.tight_layout()

    _drop_oversized_labels(fig, ax1, bar_width, texts1, values)
    _drop_oversized_labels(fig, ax2, bar_width, texts2, element_values)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor(), bbox_extra_artists=(legend1, legend2), bbox_inches="tight")
    plt.close(fig)

    return base64.b64encode(buf.getvalue()).decode("ascii")


def render_process_cost_curves(component: "Component", current_qty: int) -> str:
    """
    Plot fabrication cost vs. production quantity for every process
    overlaid on one chart. Each process gets a fixed-order categorical
    color (never cycled); within a color, dashed is an ideal material AND
    ideal geometry (Rc=1, the DFM baseline design) and solid is this
    material and this geometry's actual Rc. No title here - the page
    supplies a caption instead.

    All cost math (Pc(N), Rc, Rc * Pc) comes from
    Component.fabrication_cost_curve() - this function only plots.

    Colors are offset by 1 so each process matches its segment in the
    stacked cost chart, where palette slot 0 is reserved for "Material".
    """
    import numpy as np

    lo = max(1.0, current_qty / 100)
    hi = max(current_qty * 100, lo * 10)
    N = np.logspace(np.log10(lo), np.log10(hi), 200)

    curve_rows = component.fabrication_cost_curve(N)
    current_rows = component.fabrication_cost_curve(current_qty)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    fig.patch.set_facecolor(_SURFACE)
    ax.set_facecolor(_SURFACE)

    for i, (row, current_row) in enumerate(zip(curve_rows, current_rows)):
        proc = row["process"]
        color = _CATEGORICAL_PALETTE[i + 1] if i + 1 < len(_CATEGORICAL_PALETTE) else _OTHER_COLOR

        ax.plot(N, row["cost_ideal"], linestyle="--", linewidth=2, color=color, label=f"{proc.name} (ideal)", zorder=3)
        ax.plot(N, row["cost_selected"], linestyle="-", linewidth=2, color=color, label=proc.name, zorder=4)

        ax.plot([current_qty], [current_row["cost_selected"]], marker="o", markersize=6, color=color, zorder=5)

    ax.axvline(current_qty, color=_BASELINE, linewidth=1, linestyle=":", zorder=2)

    ax.set_xscale("log")
    ax.set_xlabel("Production quantity", color=_INK_SECONDARY, fontsize=9)
    ax.set_ylabel("Fabrication cost ($)", color=_INK_SECONDARY, fontsize=9)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"${v:,.0f}"))
    ax.tick_params(axis="both", colors=_INK_MUTED, labelsize=8)

    ax.yaxis.grid(True, color=_GRIDLINE, linewidth=1, zorder=0)
    ax.xaxis.grid(True, color=_GRIDLINE, linewidth=1, zorder=0)
    ax.set_axisbelow(True)

    for spine_name in ("top", "right"):
        ax.spines[spine_name].set_visible(False)
    ax.spines["left"].set_color(_BASELINE)
    ax.spines["bottom"].set_color(_BASELINE)

    ax.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=_INK_SECONDARY)

    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
    plt.close(fig)

    return base64.b64encode(buf.getvalue()).decode("ascii")


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

    material_fraction_rows, warnings = tea_api.material_cost_breakdown(material)

    return {
        "summary": summary,
        "chart_png_base64": render_cost_chart(summary, material_fraction_rows),
        "process_cost_curve_chart_png_base64": render_process_cost_curves(component, req.production_qty),
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
