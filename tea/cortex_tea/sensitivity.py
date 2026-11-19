"""
One-at-a-time (OAT) sensitivity analysis for the TEA cost model.

Takes the same inputs as tea_api.evaluate(), builds the base case, then
swings every model input one at a time between a low and a high value
while holding everything else at its base value, and records how the
chosen output metric moves. The result feeds plotting.tornado_png().

    from cortex_tea import sensitivity, plotting

    result = sensitivity.one_at_a_time(
        "EROFER97", ["CNC", "Hot Rolling"],
        volume_mm3=3_000_000, production_qty=100,
        rel_range=0.20,                          # default +/-20%
        overrides={"price:W": (30, 70)},         # absolute low/high
    )
    plotting.tornado_png(result, path="tornado.png")
    result.to_csv("tornado.csv")

Parameters swung (keys, as used by `overrides`, `include` and `exclude`):

    price:<el>             element price ($/kg), every element in the material
    wt:<el>                element wt%, every non-remainder element; the
                           remainder element absorbs the change
    density                material density (g/cm3)
    volume_mm3             component volume (mm3)
    production_qty         production quantity N
    Cmp/Cc/Cs/Ct/Cf:<proc> DFM relative cost coefficients, per process
    Wc:<proc>              waste coefficient, per process
    alphaT:<proc>, beta:<proc>
                           processes given alphaT/beta directly
    equipment_cost:<proc>, tooling_cost:<proc>, processing_time:<proc>
                           processes built from Low..Very High levels
    level:<proc>:tooling|equipment|time
                           +/-1 level step (only with level_steps=True)
    amortization_years, utilization, escalation_rate
                           equipment amortization constants (only when a
                           level-built process is present)

Defaults: continuous inputs swing base*(1 -/+ rel_range). DFM coefficients
(Cmp/Cc/Cs/Ct/Cf) and Wc are >= 1 by definition (1 = ideal), so their
default low end is clamped at 1.0; utilization is clamped at 1.0 on the
high end. An explicit override is used as given.

Composition: when wt:<el> moves, the remainder element (e.g. Fe) takes up
the difference, so a wt% swing measures the price difference between that
element and the remainder. A blended material has no remainder element, so
its wt rows are skipped.

Nothing here changes the model's equations: each case is evaluated with
tea1.Component.manufacturing_cost(), the same code path as evaluate().
"""

import contextlib
import copy
import csv
import fnmatch
import io
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from . import cost_variables
from . import material_def
from . import process_def
from . import tea_api
from .tea1 import Component

METRICS = {
    "Total cost": "$ per part",
    "Unit cost": "$/kg",
    "Material cost": "$ per part",
    "Processing cost": "$ per part",
}

PRICE_SOURCE = "cortex_tea/cost_variables.py material_cost_database (cost year and source not stated)"

# Base values of the amortization constants in process_def.compute_alphaT_beta:
#   alphaT = equipment_cost / (years*24*365*utilization) * (1+rate)**years * processing_time
_BASE_YEARS = 5.0
_BASE_UTILIZATION = 0.5
_BASE_ESCALATION = 0.05

_LEVEL_TABLES = {
    "tooling": process_def.TOOLING_COST_BY_LEVEL,
    "equipment": process_def.EQUIPMENT_COST_BY_LEVEL,
    "time": process_def.PROCESSING_TIME_BY_LEVEL,
}
_LEVEL_FIELD = {"tooling": "tooling_cost", "equipment": "equipment_cost", "time": "processing_time"}

_COEFF_NAMES = ("Cmp", "Cc", "Cs", "Ct", "Cf", "Wc")


class InvalidValue(ValueError):
    """A perturbed input value the model can't take (e.g. wt% sum > 100)."""


# ----------------------------------------------------------------------
# Case: a flat, deep-copyable description of one model evaluation


def _alphaT_from_levels(proc: dict, amort: dict) -> float:
    hours = amort["years"] * 24 * 365 * amort["utilization"]
    escalation = (1 + amort["escalation_rate"]) ** amort["years"]
    return proc["equipment_cost"] / hours * escalation * proc["processing_time"]


def _case_from_component(component: Component) -> dict:
    mat = component.material
    rem = [el for el, d in mat.composition.items() if d["type"] == "rem"]
    prices = {
        el: cost_variables.material_cost_database[el].cost
        for el in mat.composition
        if el in cost_variables.material_cost_database
    }
    amort = {"years": _BASE_YEARS, "utilization": _BASE_UTILIZATION, "escalation_rate": _BASE_ESCALATION}

    procs = []
    for p in component.processes:
        entry = {"name": p.name, "mode": "direct", "alphaT": p.alphaT, "beta": p.beta}
        if p.tooling_level and p.equipment_level and p.time_level:
            levels = {"tooling": p.tooling_level, "equipment": p.equipment_level, "time": p.time_level}
            candidate = dict(entry, mode="levels", levels=levels)
            for kind, level in levels.items():
                candidate[_LEVEL_FIELD[kind]] = _LEVEL_TABLES[kind][level]
            # Only treat it as level-built if the levels reproduce its alphaT/beta.
            a = _alphaT_from_levels(candidate, amort)
            if abs(a - p.alphaT) <= 1e-9 * abs(p.alphaT) and candidate["tooling_cost"] == p.beta:
                entry = candidate
        procs.append(entry)

    return {
        "material_name": mat.name,
        "density": mat.density,
        "composition": copy.deepcopy(mat.composition),
        "remainder": rem[0] if len(rem) == 1 else None,
        "prices": prices,
        "procs": procs,
        "Cmp": [mat.get_Cmp(p.name) for p in component.processes],
        "Cc": list(component.Cc),
        "Cs": list(component.Cs),
        "Ct": list(component.Ct),
        "Cf": list(component.Cf),
        "Wc": list(component.Wc),
        "amort": amort,
        "volume_mm3": component.Vf,
        "production_qty": component.N,
        "component_name": component.name,
    }


def _material_cost_per_kg(case: dict) -> float:
    """$/kg from composition, same rule as Material.calculate_cost_breakdown:
    'alloy' and 'rem' elements pay their price, 'residual' elements are free,
    elements with no price contribute 0."""
    total = 0.0
    for el, d in case["composition"].items():
        if d["type"] in ("alloy", "rem") and el in case["prices"]:
            total += d["wt"] / 100.0 * case["prices"][el]
    return total


def _component_from_case(case: dict) -> Component:
    mat = material_def.Material(case["material_name"])
    mat.composition = copy.deepcopy(case["composition"])
    mat.density = case["density"]
    mat.cost = _material_cost_per_kg(case)
    for p, cmp_value in zip(case["procs"], case["Cmp"]):
        mat.Cmp_map[p["name"]] = {"Cmp": cmp_value}

    processes = []
    for p in case["procs"]:
        if p["mode"] == "levels":
            alphaT, beta = _alphaT_from_levels(p, case["amort"]), p["tooling_cost"]
        else:
            alphaT, beta = p["alphaT"], p["beta"]
        processes.append(process_def.Process(p["name"], alphaT, beta))

    return Component(
        name=case["component_name"],
        volume_mm3=case["volume_mm3"],
        material=mat,
        processes=processes,
        production_qty=case["production_qty"],
        Cc=list(case["Cc"]), Cs=list(case["Cs"]), Ct=list(case["Ct"]),
        Cf=list(case["Cf"]), Wc=list(case["Wc"]),
    )


def _evaluate(case: dict, metric: str) -> float:
    return _component_from_case(case).manufacturing_cost()[metric]


# ----------------------------------------------------------------------
# Parameters


@dataclass
class Parameter:
    key: str
    label: str
    group: str
    units: str
    base: object
    low: object
    high: object
    apply: Callable[[dict, object], None] = field(repr=False)
    note: str = ""


def _require_positive(value, what):
    if value <= 0:
        raise InvalidValue(f"{what} must be > 0 (got {value:g})")


def _set_top(name, what):
    def apply(case, value):
        _require_positive(value, what)
        case[name] = value
    return apply


def _set_price(el):
    def apply(case, value):
        if value < 0:
            raise InvalidValue(f"price of {el} must be >= 0 (got {value:g})")
        case["prices"][el] = value
    return apply


def _set_wt(el):
    def apply(case, value):
        _require_positive(value, f"wt% of {el}")
        comp = case["composition"]
        rem = case["remainder"]
        comp[el]["wt"] = value
        others = sum(d["wt"] for e, d in comp.items() if e != rem)
        if others >= 100:
            raise InvalidValue(f"listed elements sum to {others:.3f} wt% (>= 100), no room for remainder {rem}")
        comp[rem]["wt"] = 100.0 - others
    return apply


def _set_coeff(name, i):
    def apply(case, value):
        _require_positive(value, name)
        case[name][i] = value
    return apply


def _set_proc(i, fld):
    def apply(case, value):
        _require_positive(value, fld)
        case["procs"][i][fld] = value
    return apply


def _set_level(i, kind):
    def apply(case, value):
        p = case["procs"][i]
        p["levels"][kind] = value
        p[_LEVEL_FIELD[kind]] = _LEVEL_TABLES[kind][value]
    return apply


def _set_amort(name):
    def apply(case, value):
        _require_positive(value, name)
        if name == "utilization" and value > 1:
            raise InvalidValue(f"utilization must be <= 1 (got {value:g})")
        case["amort"][name] = value
    return apply


def _set_escalation(case, value):
    if value <= -1:
        raise InvalidValue(f"escalation_rate must be > -1 (got {value:g})")
    case["amort"]["escalation_rate"] = value


def _build_parameters(case: dict, rel_range: float, level_steps: bool) -> List[Parameter]:
    r = rel_range

    def span(base, floor=None, ceil=None):
        lo, hi = base * (1 - r), base * (1 + r)
        if floor is not None:
            lo = max(lo, floor) if base >= floor else lo
        if ceil is not None:
            hi = min(hi, ceil) if base <= ceil else hi
        return lo, hi

    params: List[Parameter] = []
    rem = case["remainder"]

    # Element prices
    for el, d in case["composition"].items():
        if el not in case["prices"]:
            continue
        base = case["prices"][el]
        note = ""
        if d["type"] == "residual":
            note = "residual element: priced at $0 in the model, so its price has no effect"
        elif d["type"] == "rem":
            note = "remainder element"
        lo, hi = span(base)
        params.append(Parameter(f"price:{el}", f"{el} price", "Element price", "$/kg", base, lo, hi, _set_price(el), note))

    # Composition
    for el, d in case["composition"].items():
        if el == rem:
            continue
        base = d["wt"]
        lo, hi = span(base)
        note = f"{rem} absorbs the change" if rem else ""
        params.append(Parameter(f"wt:{el}", f"{el} wt%", "Composition", "wt%", base, lo, hi, _set_wt(el), note))

    # Material and part
    lo, hi = span(case["density"])
    mass_note = "density and volume both scale part mass only; same relative swing gives the same result"
    params.append(Parameter("density", "Density", "Material", "g/cm3", case["density"], lo, hi,
                            _set_top("density", "density"), mass_note))
    lo, hi = span(case["volume_mm3"])
    params.append(Parameter("volume_mm3", "Volume", "Part", "mm3", case["volume_mm3"], lo, hi,
                            _set_top("volume_mm3", "volume"), mass_note))
    lo, hi = span(case["production_qty"])
    params.append(Parameter("production_qty", "Production quantity", "Part", "parts", case["production_qty"], lo, hi,
                            _set_top("production_qty", "production quantity")))

    # Per-process coefficients and cost terms
    any_levels = False
    for i, p in enumerate(case["procs"]):
        name = p["name"]
        ct, cf = case["Ct"][i], case["Cf"][i]
        smaller = None if ct == cf else ("Ct" if ct < cf else "Cf")
        for coeff in _COEFF_NAMES:
            base = case[coeff][i]
            lo, hi = span(base, floor=1.0)
            if coeff == "Wc":
                note = "all Wc multiply into material mass"
            elif coeff in ("Ct", "Cf") and smaller is None:
                note = "Ct = Cf and only max(Ct, Cf) enters Rc: raising either raises Rc, lowering either has no effect"
            elif coeff == smaller:
                other = "Cf" if coeff == "Ct" else "Ct"
                note = f"only max(Ct, Cf) enters Rc; {coeff} is the smaller, no effect until it passes {other}"
            else:
                note = f"same multiplier on {name} Rc as the other Rc coefficients"
            params.append(Parameter(f"{coeff}:{name}", f"{coeff} ({name})", "DFM coefficient", "-", base, lo, hi,
                                    _set_coeff(coeff, i), note))

        if p["mode"] == "levels":
            any_levels = True
            for fld, label, units in (
                ("equipment_cost", "equipment cost", "$"),
                ("tooling_cost", "tooling cost", "$"),
                ("processing_time", "processing time", "level value"),
            ):
                base = p[fld]
                lo, hi = span(base)
                note = f"same multiplier on {name} alphaT as {'processing time' if fld == 'equipment_cost' else 'equipment cost'}" \
                    if fld != "tooling_cost" else f"= beta for {name}"
                params.append(Parameter(f"{fld}:{name}", f"{name} {label}", "Process", units, base, lo, hi,
                                        _set_proc(i, fld), note))
            if level_steps:
                for kind, level in p["levels"].items():
                    idx = process_def.COST_LEVELS.index(level)
                    lo_lvl = process_def.COST_LEVELS[max(idx - 1, 0)]
                    hi_lvl = process_def.COST_LEVELS[min(idx + 1, len(process_def.COST_LEVELS) - 1)]
                    params.append(Parameter(f"level:{name}:{kind}", f"{name} {kind} level", "Process level", "level",
                                            level, lo_lvl, hi_lvl, _set_level(i, kind),
                                            f"same input as {_LEVEL_FIELD[kind]}:{name}, stepped by one level"))
        else:
            for fld, label in (("alphaT", "alphaT"), ("beta", "beta")):
                base = p[fld]
                lo, hi = span(base)
                params.append(Parameter(f"{fld}:{name}", f"{name} {label}", "Process", "$", base, lo, hi, _set_proc(i, fld)))

    if any_levels:
        a = case["amort"]
        lo, hi = span(a["years"])
        params.append(Parameter("amortization_years", "Amortization years", "Amortization", "yr", a["years"], lo, hi,
                                _set_amort("years"), "also the escalation period"))
        lo, hi = span(a["utilization"], ceil=1.0)
        params.append(Parameter("utilization", "Equipment utilization", "Amortization", "fraction", a["utilization"], lo, hi,
                                _set_amort("utilization")))
        lo, hi = span(a["escalation_rate"])
        params.append(Parameter("escalation_rate", "Escalation rate", "Amortization", "fraction/yr", a["escalation_rate"],
                                lo, hi, _set_escalation))

    return params


def _filter(params, include, exclude):
    def matches(key, patterns):
        return any(fnmatch.fnmatchcase(key, pat) for pat in patterns)

    if include:
        params = [p for p in params if matches(p.key, include)]
    if exclude:
        params = [p for p in params if not matches(p.key, exclude)]
    return params


# ----------------------------------------------------------------------
# Result


@dataclass
class SensitivityResult:
    metric: str
    units: str
    base_value: float
    rows: List[dict]
    basis: dict

    def to_csv(self, path) -> None:
        cols = ["key", "label", "group", "units", "base_input", "low_input", "high_input",
                "low_output", "high_output", "low_delta", "high_delta", "swing", "note"]
        with open(path, "w", newline="", encoding="utf-8") as fh:
            for k, v in self.basis.items():
                fh.write(f"# {k}: {v}\n")
            writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.rows)

    def to_dict(self) -> dict:
        return {"metric": self.metric, "units": self.units, "base_value": self.base_value,
                "rows": self.rows, "basis": self.basis}


def _model_version() -> str:
    try:
        from importlib.metadata import version
        return f"cortex-tea {version('cortex-tea')}"
    except Exception:
        return "cortex-tea (version unknown)"


def _base_case(material, processes, volume_mm3, production_qty, geometry, material_dir, process_dir, geometry_dir):
    if isinstance(material, material_def.Material):
        material = copy.deepcopy(material)  # build_component may add Cmp entries
    with contextlib.redirect_stdout(io.StringIO()):
        component = tea_api.make_component(
            material, processes, volume_mm3, production_qty,
            geometry=geometry, material_dir=material_dir,
            process_dir=process_dir, geometry_dir=geometry_dir,
        )
    return component, _case_from_component(component)


def list_parameters(material, processes, volume_mm3, production_qty, *, geometry=None,
                    rel_range: float = 0.20, level_steps: bool = False,
                    material_dir=None, process_dir=None, geometry_dir=None) -> List[dict]:
    """Every parameter one_at_a_time() would swing, with its base value and
    default low/high, without evaluating anything."""
    _, case = _base_case(material, processes, volume_mm3, production_qty, geometry,
                         material_dir, process_dir, geometry_dir)
    return [
        {"key": p.key, "label": p.label, "group": p.group, "units": p.units,
         "base": p.base, "low": p.low, "high": p.high, "note": p.note}
        for p in _build_parameters(case, rel_range, level_steps)
    ]


def one_at_a_time(material, processes, volume_mm3, production_qty, *, geometry=None,
                  metric: str = "Total cost", rel_range: float = 0.20,
                  overrides: Optional[Dict[str, Tuple[object, object]]] = None,
                  include: Optional[List[str]] = None, exclude: Optional[List[str]] = None,
                  level_steps: bool = False, component_name: str = "Component",
                  material_dir=None, process_dir=None, geometry_dir=None) -> SensitivityResult:
    """Swing each input between its low and high value, one at a time.

    Same material/process/geometry inputs as tea_api.evaluate(). Returns a
    SensitivityResult whose rows are sorted by swing, largest first. Each
    row: key, label, group, units, base_input, low_input, high_input,
    low_output, high_output, low_delta, high_delta, swing
    (= |high_output - low_output|, a skipped side counted at the base value),
    note. A side whose input the model can't take is recorded as
    low_output/high_output = None with the reason in note.
    """
    component, _ = _base_case(material, processes, volume_mm3, production_qty, geometry,
                              material_dir, process_dir, geometry_dir)
    component.name = component_name
    return from_component(component, metric=metric, rel_range=rel_range, overrides=overrides,
                          include=include, exclude=exclude, level_steps=level_steps)


def from_component(component: Component, *, metric: str = "Total cost", rel_range: float = 0.20,
                   overrides: Optional[Dict[str, Tuple[object, object]]] = None,
                   include: Optional[List[str]] = None, exclude: Optional[List[str]] = None,
                   level_steps: bool = False) -> SensitivityResult:
    """one_at_a_time() for an already-built Component (e.g. from
    tea_api.build_component, as the web GUI does). The component is not
    modified."""
    if metric not in METRICS:
        raise tea_api.TeaError(f"metric must be one of {list(METRICS)}")
    if not 0 < rel_range < 1:
        raise tea_api.TeaError("rel_range must be between 0 and 1 (e.g. 0.2 for +/-20%)")

    case = _case_from_component(component)

    base_value = _evaluate(case, metric)
    reference = component.manufacturing_cost()[metric]
    if abs(base_value - reference) > 1e-9 * max(1.0, abs(reference)):
        raise RuntimeError(
            f"Sensitivity base case ({base_value}) does not reproduce the model ({reference}); "
            "the case reconstruction is out of step with tea1.Component."
        )

    params = _build_parameters(case, rel_range, level_steps)
    overrides = dict(overrides or {})
    known = {p.key for p in params}
    unknown = [k for k in overrides if k not in known]
    if unknown:
        raise tea_api.TeaError(f"Unknown override key(s) {unknown}. Use list_parameters() to see valid keys.")
    for p in params:
        if p.key in overrides:
            p.low, p.high = overrides[p.key]
    params = _filter(params, include, exclude)

    rows = []
    for p in params:
        outputs, notes = {}, [p.note] if p.note else []
        for side, value in (("low", p.low), ("high", p.high)):
            trial = copy.deepcopy(case)
            try:
                p.apply(trial, value)
                outputs[side] = _evaluate(trial, metric)
            except (InvalidValue, KeyError, TypeError) as exc:
                outputs[side] = None
                notes.append(f"{side} skipped: {exc}")
        lo_out = outputs["low"] if outputs["low"] is not None else base_value
        hi_out = outputs["high"] if outputs["high"] is not None else base_value
        rows.append({
            "key": p.key, "label": p.label, "group": p.group, "units": p.units,
            "base_input": p.base, "low_input": p.low, "high_input": p.high,
            "low_output": outputs["low"], "high_output": outputs["high"],
            "low_delta": None if outputs["low"] is None else outputs["low"] - base_value,
            "high_delta": None if outputs["high"] is None else outputs["high"] - base_value,
            "swing": abs(hi_out - lo_out),
            "note": "; ".join(n for n in notes if n),
        })
    rows.sort(key=lambda row: row["swing"], reverse=True)

    basis = {
        "model": _model_version(),
        "material": case["material_name"],
        "processes": ", ".join(p["name"] for p in case["procs"]),
        "geometry": getattr(component.geometry, "name", None) or "none (coefficients from overrides/defaults)",
        "volume_mm3": case["volume_mm3"],
        "production_qty": case["production_qty"],
        "metric": f"{metric} ({METRICS[metric]})",
        "base_value": base_value,
        "default_range": f"+/-{rel_range:.0%} (DFM coefficients and Wc floored at 1.0)"
                         + ("; level rows are +/-1 level step, a much larger change" if level_steps else ""),
        "overrides": ", ".join(f"{k}={v[0]}..{v[1]}" for k, v in overrides.items()) or "none",
        "element_prices": PRICE_SOURCE,
        "method": "one at a time; all other inputs at base",
    }
    return SensitivityResult(metric, METRICS[metric], base_value, rows, basis)
