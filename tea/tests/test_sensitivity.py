"""Tests for cortex_tea.sensitivity. Each test checks a property with a
known answer from the model equations, not just that the code runs."""

import math

import pytest

from cortex_tea import cost_variables, plotting, sensitivity, tea_api

MAT = "EROFER97"
PROCS = ["CNC", "Hot Rolling"]
VOL = 3_000_000
N = 100


def _run(**kw):
    args = dict(material=MAT, processes=PROCS, volume_mm3=VOL, production_qty=N)
    args.update(kw)
    return sensitivity.one_at_a_time(**args)


def _row(result, key):
    return next(r for r in result.rows if r["key"] == key)


def _summary():
    return tea_api.evaluate(MAT, PROCS, volume_mm3=VOL, production_qty=N)["summary"]


@pytest.mark.parametrize("metric", list(sensitivity.METRICS))
def test_base_value_matches_evaluate(metric):
    assert _run(metric=metric).base_value == pytest.approx(_summary()[metric], rel=1e-12)


def test_cmp_swing_is_linear_in_process_cost():
    summary = _summary()
    cnc_cost = next(p["Cost"] for p in summary["Processes"] if p["Process"] == "CNC")
    row = _row(_run(), "Cmp:CNC")
    assert row["base_input"] == 4
    assert row["swing"] == pytest.approx(0.4 * cnc_cost, rel=1e-12)


def test_production_qty_matches_beta_over_n():
    result = _run()
    row = _row(result, "production_qty")
    component = tea_api.make_component(MAT, PROCS, volume_mm3=VOL, production_qty=N)
    expected = sum(
        component.get_Rc(i) * p.beta * (1 / row["low_input"] - 1 / row["high_input"])
        for i, p in enumerate(component.processes)
    )
    assert row["low_output"] - row["high_output"] == pytest.approx(expected, rel=1e-12)


def test_residual_element_price_has_no_effect():
    result = _run()
    for el in ("Ni", "Si", "C"):  # residuals in EUROFER97
        assert _row(result, f"price:{el}")["swing"] == 0


def test_wt_change_is_absorbed_by_remainder():
    result = _run(metric="Material cost")
    row = _row(result, "wt:Cr")
    delta_wt = row["high_input"] - row["base_input"]
    db = cost_variables.material_cost_database
    mass = 8 * VOL / 1e6
    expected = delta_wt / 100 * (db["Cr"].cost - db["Fe"].cost) * mass * 1.0  # Wc = 1
    assert row["high_delta"] == pytest.approx(expected, rel=1e-12)

    # The composition stays at 100 wt% in the perturbed case.
    _, case = sensitivity._base_case(MAT, PROCS, VOL, N, None, None, None, None)
    trial = sensitivity.copy.deepcopy(case)
    sensitivity._set_wt("Cr")(trial, row["high_input"])
    assert sum(d["wt"] for d in trial["composition"].values()) == pytest.approx(100.0, abs=1e-12)
    assert trial["composition"]["Fe"]["wt"] == pytest.approx(case["composition"]["Fe"]["wt"] - delta_wt, abs=1e-12)


def test_database_material_not_mutated():
    before = _summary()["Total cost"]
    _run()
    _run(overrides={"Cmp:CNC": (1.0, 10.0)})
    assert _summary()["Total cost"] == before


def test_smaller_of_ct_cf_has_no_effect():
    procs = [{"name": "CNC", "Ct": 1.0, "Cf": 2.0}, "Hot Rolling"]
    result = _run(processes=procs)
    assert _row(result, "Ct:CNC")["swing"] == 0
    assert _row(result, "Cf:CNC")["swing"] > 0


def test_coefficient_low_end_floored_at_one():
    result = _run()
    row = _row(result, "Cc:CNC")  # base 1.0
    assert row["low_input"] == 1.0
    assert row["low_delta"] == 0


def test_invalid_side_is_skipped_not_raised():
    result = _run(overrides={"wt:Cr": (9.0, 99.0)})
    row = _row(result, "wt:Cr")
    assert row["high_output"] is None
    assert "high skipped" in row["note"]


def test_unknown_override_key_raises():
    with pytest.raises(tea_api.TeaError):
        _run(overrides={"price:Unobtainium": (1, 2)})


def test_level_built_processes_expose_amortization():
    keys = {p["key"] for p in sensitivity.list_parameters(MAT, PROCS, VOL, N)}
    assert {"equipment_cost:CNC", "tooling_cost:CNC", "processing_time:CNC", "amortization_years",
            "utilization", "escalation_rate"} <= keys
    assert "alphaT:CNC" not in keys


def test_equipment_cost_and_time_are_same_multiplier():
    result = _run()
    assert _row(result, "equipment_cost:CNC")["swing"] == pytest.approx(
        _row(result, "processing_time:CNC")["swing"], rel=1e-12)


def test_level_steps():
    result = _run(level_steps=True)
    row = _row(result, "level:CNC:time")
    assert (row["low_input"], row["base_input"], row["high_input"]) == ("Medium", "Medium-High", "High")
    assert row["high_output"] > result.base_value > row["low_output"]


def test_include_exclude():
    result = _run(include=["price:*", "wt:*"], exclude=["wt:*"])
    assert result.rows and all(r["key"].startswith("price:") for r in result.rows)


def test_rows_sorted_and_csv(tmp_path):
    result = _run()
    swings = [r["swing"] for r in result.rows]
    assert swings == sorted(swings, reverse=True)
    path = tmp_path / "t.csv"
    result.to_csv(path)
    text = path.read_text(encoding="utf-8")
    assert "# element_prices:" in text and "Cmp:CNC" in text


def test_tornado_png(tmp_path):
    png = plotting.tornado_png(_run(), path=tmp_path / "t.png")
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert (tmp_path / "t.png").read_bytes() == png
    assert math.isfinite(len(png))
