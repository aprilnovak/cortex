"""
Programmatic Python API for the TEA cost model.

Lets a researcher import this module directly (no running server, no
interactive prompts) and run the cost model on one material or screen a
batch of candidate materials, getting back plain JSON-serializable dicts.

    import sys
    sys.path.insert(0, "/path/to/cortex/tea")
    import tea_api

    result = tea_api.evaluate("EROFER97", ["CNC", "Hot Rolling"],
                               volume_mm3=3_000_000, production_qty=100)

    candidates = ["EROFER97", {"name": "V1", "density": 7.6,
                  "composition": {"Cr": {"wt": 9.0, "type": "alloy"}},
                  "remainder_element": "Fe"}]
    results = tea_api.screen(candidates, ["CNC"],
                              volume_mm3=1_000_000, production_qty=50)

    # Or run every material file in a directory (e.g. a folder of candidate
    # alloy .py files, same format as materials_database/):
    results = tea_api.screen_directory("candidate_alloys/", ["CNC"],
                                        volume_mm3=1_000_000, production_qty=50)

This module is also the shared backend for tea/webapp/server.py, which
wraps the TeaError/TeaNotFoundError exceptions raised here into HTTP
responses for the browser GUI.
"""

import contextlib
import importlib.util
import io
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Union

import material_def
import process_def
import component_def
import cost_variables
from tea1 import Component

_TEA_DIR = Path(__file__).resolve().parent
MATERIALS_DB_DIR = str(_TEA_DIR / "materials_database")
PROCESSES_DB_DIR = str(_TEA_DIR / "processes_database")
GEOMETRIES_DB_DIR = str(_TEA_DIR / "geometries_database")


class TeaError(ValueError):
    """Invalid input: a material/process/geometry failed to build, or a
    required value (e.g. Cmp) was missing and had no fallback."""


class TeaNotFoundError(TeaError):
    """A named material/process/geometry wasn't found in the database or
    the supplied extra directory."""


# ----------------------------------------------------------------------
# Helpers for translating material_def/process_def/component_def's
# "return None + print an ASCII error box to stdout" convention into a
# raised exception with a plain-text message.


def _captured_call(fn, *args, **kwargs):
    """Call fn, capturing anything it prints to stdout."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = fn(*args, **kwargs)
    return result, buf.getvalue()


def _extract_error(captured_text: str, fallback: str) -> str:
    """Pull the message out of print_error's ASCII error box, if present."""
    lines = []
    in_box = False
    for line in captured_text.splitlines():
        if "ERROR" in line and set(line.strip()) & {"-"}:
            in_box = True
            continue
        if not in_box:
            continue
        stripped = line.strip()
        if stripped and set(stripped) <= {"-"}:
            break
        if stripped.startswith("|") and stripped.endswith("|"):
            content = stripped.strip("|").strip()
            if content:
                lines.append(content)
    text = re.sub(r"\s+", " ", " ".join(lines)).strip()
    return text or fallback


# ----------------------------------------------------------------------
# Database + optional extra-directory lookups


def _load_material_lookup(material_dir: Optional[str] = None) -> Dict[str, "material_def.Material"]:
    lookup = {}
    db_materials, _ = _captured_call(material_def.load_materials, MATERIALS_DB_DIR)
    lookup.update({m.name: m for m in db_materials if m is not None})

    if material_dir:
        if not os.path.isdir(material_dir):
            raise TeaError(f"Material directory '{material_dir}' not found.")
        custom, _ = _captured_call(material_def.load_materials, material_dir)
        lookup.update({m.name: m for m in custom if m is not None})

    return lookup


def _load_process_lookup(process_dir: Optional[str] = None) -> Dict[str, "process_def.Process"]:
    lookup = {}
    try:
        db_processes, _ = _captured_call(process_def.load_processes, PROCESSES_DB_DIR)
    except FileNotFoundError:
        db_processes = []
    lookup.update({p.name: p for p in db_processes if p is not None})

    if process_dir:
        if not os.path.isdir(process_dir):
            raise TeaError(f"Process directory '{process_dir}' not found.")
        custom, _ = _captured_call(process_def.load_processes, process_dir)
        lookup.update({p.name: p for p in custom if p is not None})

    return lookup


def _load_geometry_lookup(geometry_dir: Optional[str] = None) -> Dict[str, "component_def.Geometry"]:
    lookup = {}
    db_geometries, _ = _captured_call(component_def.load_geometries, GEOMETRIES_DB_DIR)
    lookup.update({g.name: g for g in db_geometries if g is not None})

    if geometry_dir:
        if not os.path.isdir(geometry_dir):
            raise TeaError(f"Geometry directory '{geometry_dir}' not found.")
        custom, _ = _captured_call(component_def.load_geometries, geometry_dir)
        lookup.update({g.name: g for g in custom if g is not None})

    return lookup


# ----------------------------------------------------------------------
# Serialization (JSON-safe dicts for a Material/Process/Geometry)


def serialize_material(m: "material_def.Material") -> dict:
    return {
        "name": m.name,
        "density": m.density,
        "composition": m.composition,
        "cmp_map": {proc: v["Cmp"] for proc, v in m.Cmp_map.items()},
        "cost": m.cost,
    }


def serialize_process(p: "process_def.Process") -> dict:
    return {
        "name": p.name,
        "alphaT": p.alphaT,
        "beta": p.beta,
        "description": p.description,
        "tooling_level": p.tooling_level,
        "equipment_level": p.equipment_level,
        "time_level": p.time_level,
    }


def serialize_geometry(g: "component_def.Geometry") -> dict:
    return {
        "name": g.name,
        "volume_mm3": g.volume_mm3,
        "section_thickness_mm": g.section_thickness_mm,
        "shape_class": g.shape_class,
        "tolerance_mm": g.tolerance_mm,
        "surface_finish_um_ra": g.surface_finish_um_ra,
        "Cc_map": g.Cc_map,
        "Cs_map": g.Cs_map,
        "Ct_map": g.Ct_map,
        "Cf_map": g.Cf_map,
    }


def material_cost_breakdown(material: "material_def.Material") -> tuple:
    """Each element's absolute ($/kg) and relative (share of the material's
    total $/kg cost) contribution, plus any warnings (e.g. missing cost
    data)."""
    breakdown, missing_elements = material.calculate_cost_breakdown(cost_variables.material_cost_database)
    rows = [
        {
            "element": row["element"],
            "wt": row["wt"],
            "type": row["type"],
            "dollar_per_kg": row["dollar_per_kg"],
            "fraction": row["dollar_per_kg"] / material.cost if material.cost else 0.0,
        }
        for row in breakdown
    ]
    rows.sort(key=lambda r: r["fraction"], reverse=True)

    warnings = []
    if missing_elements:
        warnings.append(f"Missing cost data for elements {missing_elements} in '{material.name}'.")
    return rows, warnings


# ----------------------------------------------------------------------
# Resolution: turn a plain str/dict/object spec into a live
# Material/Process/Geometry.

MaterialSpec = Union[str, dict, "material_def.Material"]
ProcessSpec = Union[str, dict, "process_def.Process"]
GeometrySpec = Union[str, "component_def.Geometry", None]


def resolve_material(material: MaterialSpec, material_dir: Optional[str] = None) -> "material_def.Material":
    """Resolve a material spec to a Material instance.

    material may be:
      - a name string (looked up in the database + optional material_dir)
      - a dict with 'composition' (custom material, via build_material):
        {"name", "density", "composition": {el: {"wt", "type"}},
         "remainder_element", "cmp_map"?}
      - a dict with 'material_a'/'material_b'/'fraction_a' (blend, via
        blend_materials); material_a/material_b may themselves be a name
        string or a custom-composition dict, resolved recursively.
      - an already-built Material instance (passthrough).
    """
    if isinstance(material, material_def.Material):
        return material

    if isinstance(material, str):
        lookup = _load_material_lookup(material_dir)
        found = lookup.get(material)
        if found is None:
            raise TeaNotFoundError(f"Material '{material}' not found.")
        return found

    if isinstance(material, dict):
        if "composition" in material:
            name = material.get("name")
            density = material.get("density")
            composition = material.get("composition")
            remainder_element = material.get("remainder_element")
            cmp_map = material.get("cmp_map", {})
            if density is None or not composition or not remainder_element:
                raise TeaError("Custom material requires density, composition, and a remainder element.")

            built, captured = _captured_call(
                material_def.build_material,
                name=name,
                density=density,
                composition=composition,
                remainder_element=remainder_element,
                Cmp_map=cmp_map,
            )
            if built is None:
                raise TeaError(_extract_error(captured, f"Failed to build material '{name}'."))
            return built

        if "material_a" in material or "material_b" in material:
            material_a_spec = material.get("material_a")
            material_b_spec = material.get("material_b")
            fraction_a = material.get("fraction_a")
            basis = material.get("basis", "volume")
            name = material.get("name")
            if material_a_spec is None or material_b_spec is None or fraction_a is None:
                raise TeaError("Blended material requires material_a, material_b, and fraction_a.")

            material_a = resolve_material(material_a_spec, material_dir)
            material_b = resolve_material(material_b_spec, material_dir)
            blended, captured = _captured_call(
                material_def.blend_materials,
                name=name or f"{material_a.name}/{material_b.name} Blend",
                material_a=material_a,
                material_b=material_b,
                fraction_a=fraction_a,
                basis=basis,
            )
            if blended is None:
                raise TeaError(_extract_error(captured, "Failed to blend materials."))
            return blended

        raise TeaError(
            "Material dict must contain 'composition' (custom material) or "
            "'material_a'/'material_b'/'fraction_a' (blend)."
        )

    raise TeaError(f"Unsupported material spec: {material!r}")


def resolve_process(process: ProcessSpec, process_dir: Optional[str] = None) -> "process_def.Process":
    """Resolve a process spec to a Process instance.

    process may be:
      - a name string (looked up in the database + optional process_dir)
      - a dict {"name", "tooling_level"/"equipment_level"/"time_level"} or
        {"name", "alphaT", "beta"}, plus optional "description" (custom
        process, via build_process)
      - an already-built Process instance (passthrough).

    Per-process coefficient overrides (Cc/Cs/Ct/Cf/Wc/Cmp) are read by the
    caller (build_component/evaluate/screen), not resolved here.
    """
    if isinstance(process, process_def.Process):
        return process

    if isinstance(process, str):
        lookup = _load_process_lookup(process_dir)
        found = lookup.get(process)
        if found is None:
            raise TeaNotFoundError(f"Process '{process}' not found.")
        return found

    if isinstance(process, dict):
        name = process.get("name")
        if "alphaT" in process or "beta" in process:
            built, captured = _captured_call(
                process_def.build_process,
                name=name,
                alphaT=process.get("alphaT"),
                beta=process.get("beta"),
                description=process.get("description", ""),
                verbose=False,
            )
        elif "tooling_level" in process or "equipment_level" in process or "time_level" in process:
            built, captured = _captured_call(
                process_def.build_process,
                name=name,
                tooling_level=process.get("tooling_level"),
                equipment_level=process.get("equipment_level"),
                time_level=process.get("time_level"),
                description=process.get("description", ""),
                verbose=False,
            )
        else:
            lookup = _load_process_lookup(process_dir)
            found = lookup.get(name)
            if found is None:
                raise TeaNotFoundError(f"Process '{name}' not found.")
            return found

        if built is None:
            raise TeaError(_extract_error(captured, f"Failed to build process '{name}'."))
        return built

    raise TeaError(f"Unsupported process spec: {process!r}")


def resolve_processes(processes: List[ProcessSpec], process_dir: Optional[str] = None) -> List["process_def.Process"]:
    return [resolve_process(p, process_dir) for p in processes]


def resolve_geometry(geometry: GeometrySpec, geometry_dir: Optional[str] = None) -> Optional["component_def.Geometry"]:
    if geometry is None or isinstance(geometry, component_def.Geometry):
        return geometry

    if isinstance(geometry, str):
        lookup = _load_geometry_lookup(geometry_dir)
        found = lookup.get(geometry)
        if found is None:
            raise TeaNotFoundError(f"Geometry '{geometry}' not found.")
        return found

    raise TeaError(f"Unsupported geometry spec: {geometry!r}")


def _process_overrides(process: ProcessSpec) -> dict:
    """Per-process Cc/Cs/Ct/Cf/Wc/Cmp overrides from a spec dict, defaulting
    Cc/Cs/Ct/Cf/Wc to 1.0 and Cmp to None (no manual fallback) - a name
    string or an already-built Process instance carries no overrides."""
    defaults = {"Cc": 1.0, "Cs": 1.0, "Ct": 1.0, "Cf": 1.0, "Wc": 1.0, "Cmp": None}
    if isinstance(process, dict):
        return {key: process.get(key, default) for key, default in defaults.items()}
    return defaults


def build_component(
    name: str,
    volume_mm3: float,
    production_qty: int,
    material: "material_def.Material",
    processes: List["process_def.Process"],
    process_specs: Optional[List[ProcessSpec]] = None,
    geometry: Optional["component_def.Geometry"] = None,
) -> "Component":
    """Build a Component from already-resolved Material/Process/Geometry
    objects. process_specs (the original, possibly-dict process specs
    passed to evaluate()/screen(), same length/order as processes) supplies
    any per-process Cc/Cs/Ct/Cf/Wc/Cmp overrides; a geometry's per-process
    coefficient map wins over an explicit override, which wins over the
    default of 1.0. If the material has no Cmp for a process and no
    override is given, raises TeaError."""
    if process_specs is None:
        process_specs = [None] * len(processes)

    def _resolve_coeff(coeff_name, coeff_map_attr):
        values = []
        for spec, proc in zip(process_specs, processes):
            coeff_map = getattr(geometry, coeff_map_attr) if geometry is not None else {}
            override = _process_overrides(spec)[coeff_name]
            values.append(coeff_map[proc.name] if proc.name in coeff_map else override)
        return values

    for spec, proc in zip(process_specs, processes):
        if proc.name not in material.Cmp_map:
            cmp_override = _process_overrides(spec)["Cmp"]
            if cmp_override is None:
                raise TeaError(
                    f"Material '{material.name}' has no compatibility cost (Cmp) for process "
                    f"'{proc.name}'. Please provide one via a process spec's 'Cmp' key."
                )
            material.add_Cmp(proc.name, cmp_override)

    Cc_list = _resolve_coeff("Cc", "Cc_map")
    Cs_list = _resolve_coeff("Cs", "Cs_map")
    Ct_list = _resolve_coeff("Ct", "Ct_map")
    Cf_list = _resolve_coeff("Cf", "Cf_map")
    Wc_list = [_process_overrides(spec)["Wc"] for spec in process_specs]

    return Component(
        name=name,
        volume_mm3=volume_mm3,
        material=material,
        processes=processes,
        production_qty=production_qty,
        geometry=geometry,
        Cc=Cc_list,
        Cs=Cs_list,
        Ct=Ct_list,
        Cf=Cf_list,
        Wc=Wc_list,
    )


# ----------------------------------------------------------------------
# Public single-material / batch entry points


def _evaluate_resolved(
    material: "material_def.Material",
    processes: List["process_def.Process"],
    process_specs: List[ProcessSpec],
    volume_mm3: float,
    production_qty: int,
    component_name: str,
    geometry: Optional["component_def.Geometry"],
) -> dict:
    """Shared by evaluate()/screen()/screen_directory(): build a Component
    from already-resolved objects and return the plain result dict."""
    component = build_component(
        component_name, volume_mm3, production_qty,
        material, processes, process_specs=process_specs, geometry=geometry,
    )
    summary = component.manufacturing_cost()
    composition_rows, warnings = material_cost_breakdown(material)
    return {
        "summary": summary,
        "material_composition": composition_rows,
        "warnings": warnings,
    }


def evaluate(
    material: MaterialSpec,
    processes: List[ProcessSpec],
    volume_mm3: float,
    production_qty: int,
    *,
    geometry: GeometrySpec = None,
    component_name: str = "Component",
    material_dir: Optional[str] = None,
    process_dir: Optional[str] = None,
    geometry_dir: Optional[str] = None,
) -> dict:
    """Run the cost model for one material against a set of processes.

    Returns a plain, JSON-serializable dict:
        {"summary": {...Component.manufacturing_cost() output...},
         "material_composition": [...per-element cost share rows...],
         "warnings": [...]}
    """
    if not processes:
        raise TeaError("At least one process must be specified.")

    resolved_material = resolve_material(material, material_dir)
    resolved_geometry = resolve_geometry(geometry, geometry_dir)
    resolved_processes = [resolve_process(p, process_dir) for p in processes]

    return _evaluate_resolved(
        resolved_material, resolved_processes, processes,
        volume_mm3, production_qty, component_name, resolved_geometry,
    )


def _candidate_label(candidate: MaterialSpec) -> str:
    if isinstance(candidate, str):
        return candidate
    if isinstance(candidate, dict):
        return candidate.get("name", repr(candidate))
    if isinstance(candidate, material_def.Material):
        return candidate.name
    return repr(candidate)


def screen(
    candidates: List[MaterialSpec],
    processes: List[ProcessSpec],
    volume_mm3: float,
    production_qty: int,
    *,
    geometry: GeometrySpec = None,
    component_name: str = "Component",
    material_dir: Optional[str] = None,
    process_dir: Optional[str] = None,
    geometry_dir: Optional[str] = None,
) -> List[dict]:
    """Run the cost model over a batch of candidate materials (e.g. for
    screening candidate alloy compositions), sharing the same processes/
    geometry/volume/production_qty across all of them.

    processes and geometry are resolved once up front. Each candidate is
    evaluated independently; a candidate that fails to resolve or build
    (bad composition, missing Cmp, etc.) is skipped and recorded with its
    error rather than aborting the whole batch.

    Returns a list of dicts, one per candidate, in input order:
        {"material": <label>, "ok": bool,
         "result": <evaluate() output or None>, "error": <str or None>}
    """
    if not processes:
        raise TeaError("At least one process must be specified.")

    resolved_processes = [resolve_process(p, process_dir) for p in processes]
    resolved_geometry = resolve_geometry(geometry, geometry_dir)

    results = []
    for candidate in candidates:
        label = _candidate_label(candidate)
        try:
            resolved_material = resolve_material(candidate, material_dir)
            result = _evaluate_resolved(
                resolved_material, resolved_processes, processes,
                volume_mm3, production_qty, component_name, resolved_geometry,
            )
            results.append({"material": label, "ok": True, "result": result, "error": None})
        except Exception as exc:
            results.append({"material": label, "ok": False, "result": None, "error": str(exc)})

    return results


def _load_materials_from_dir(material_dir: str):
    """Yield (name, Material) for every successfully-built material file in
    material_dir, or (name, TeaError) for one that failed to build - same
    module-attribute convention as material_def.load_materials (name,
    density, composition, remainder_element, Cmp_map), but building each
    file individually so a failure can be attributed to its specific file/
    name (material_def.load_materials only prints per-file errors, it
    doesn't return which file a None entry came from)."""
    if not os.path.isdir(material_dir):
        raise TeaError(f"Material directory '{material_dir}' not found.")

    for filename in sorted(os.listdir(material_dir)):
        if not filename.endswith(".py"):
            continue

        filepath = os.path.join(material_dir, filename)
        module_name = os.path.splitext(filename)[0]
        spec = importlib.util.spec_from_file_location(module_name, filepath)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        name = module.name
        built, captured = _captured_call(
            material_def.build_material,
            name=name,
            density=module.density,
            composition=module.composition,
            remainder_element=module.remainder_element,
            Cmp_map=getattr(module, "Cmp_map", {}),
        )
        if built is None:
            yield name, TeaError(_extract_error(captured, f"Failed to build material '{name}' from '{filename}'."))
        else:
            yield name, built


def screen_directory(
    material_dir: str,
    processes: List[ProcessSpec],
    volume_mm3: float,
    production_qty: int,
    *,
    geometry: GeometrySpec = None,
    component_name: str = "Component",
    process_dir: Optional[str] = None,
    geometry_dir: Optional[str] = None,
) -> List[dict]:
    """Run the cost model for every material file in material_dir (same
    module-attribute convention as materials_database/ - each file defines
    name/density/composition/remainder_element/Cmp_map at module level),
    sharing the same processes/geometry/volume/production_qty across all of
    them. Same skip-and-record behavior as screen(): a file that fails to
    build, or a material that fails to build a Component (e.g. missing
    Cmp), is recorded with its error rather than aborting the run.

    Returns a list of dicts in filename order, same shape as screen():
        {"material": <name>, "ok": bool,
         "result": <evaluate() output or None>, "error": <str or None>}
    """
    if not processes:
        raise TeaError("At least one process must be specified.")

    resolved_processes = [resolve_process(p, process_dir) for p in processes]
    resolved_geometry = resolve_geometry(geometry, geometry_dir)

    results = []
    for label, material_or_error in _load_materials_from_dir(material_dir):
        if isinstance(material_or_error, Exception):
            results.append({"material": label, "ok": False, "result": None, "error": str(material_or_error)})
            continue
        try:
            result = _evaluate_resolved(
                material_or_error, resolved_processes, processes,
                volume_mm3, production_qty, component_name, resolved_geometry,
            )
            results.append({"material": label, "ok": True, "result": result, "error": None})
        except Exception as exc:
            results.append({"material": label, "ok": False, "result": None, "error": str(exc)})

    return results
