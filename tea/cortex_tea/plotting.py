"""
PNG charts for TEA results, shared by the web GUI, Sirepo and scripts.

    from cortex_tea import tea_api, plotting

    component = tea_api.make_component("EROFER97", ["CNC", "Hot Rolling"],
                                       volume_mm3=3_000_000, production_qty=100)
    png = plotting.cost_breakdown_png(component, path="breakdown.png")
    png = plotting.process_cost_curves_png(component, path="curves.png")

Each function returns the PNG as bytes, and also writes it to `path` if one
is given. Figures are drawn on matplotlib's object API rather than pyplot,
so importing this module never changes the caller's matplotlib backend.
"""

import io
from pathlib import Path
from typing import Optional, Union

import matplotlib.ticker as mticker
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import ConnectionPatch

from .tea_api import material_cost_breakdown


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


def _new_figure(figsize, dpi) -> Figure:
    """A Figure with its own Agg canvas, independent of pyplot's state."""
    fig = Figure(figsize=figsize, dpi=dpi)
    FigureCanvasAgg(fig)
    return fig


def _save_png(fig: Figure, path, **savefig_kwargs) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", **savefig_kwargs)
    png = buf.getvalue()
    if path is not None:
        Path(path).write_bytes(png)
    return png


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


def cost_breakdown_png(component, path: Optional[Union[str, Path]] = None) -> bytes:
    """
    Cost breakdown chart for a Component (see tea_api.make_component).

    Stacked bar of Material + each process's cost, with a "zoom" panel next
    to it showing the Material segment's cost broken down by element -
    shaded in the same blue as the Material segment (light -> dark, a
    sequential one-hue ramp) and linked to it with connector lines, so it
    reads as a detail view of that one segment (matplotlib's "bar of pie"
    pattern, bar-to-bar instead of pie-to-bar).
    """
    summary = component.manufacturing_cost()
    material_fraction_rows, _ = material_cost_breakdown(component.material)

    labels = ["Material"] + [p["Process"] for p in summary["Processes"]]
    values = [summary["Cost Breakdown"][label] for label in labels]
    colors = [_CATEGORICAL_PALETTE[i] if i < len(_CATEGORICAL_PALETTE) else _OTHER_COLOR for i in range(len(labels))]

    visible = [r for r in material_fraction_rows if r["fraction"] > 0]
    element_labels = [r["element"] for r in visible]
    element_values = [r["fraction"] * 100 for r in visible]
    element_colors = [_BLUE_RAMP[i] if i < len(_BLUE_RAMP) else _OTHER_COLOR for i in range(len(element_labels))]

    bar_width = 1.2
    fig = _new_figure(figsize=(8.6, 5.4), dpi=150)
    ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.3, 1]})
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

    return _save_png(
        fig, path, facecolor=fig.get_facecolor(), bbox_extra_artists=(legend1, legend2), bbox_inches="tight",
    )


def process_cost_curves_png(
    component, production_qty: Optional[int] = None, path: Optional[Union[str, Path]] = None,
) -> bytes:
    """
    Fabrication cost vs. production quantity for a Component (see
    tea_api.make_component). The current-quantity marker defaults to the
    component's own production quantity.

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
    current_qty = component.N if production_qty is None else production_qty

    lo = max(1.0, current_qty / 100)
    hi = max(current_qty * 100, lo * 10)
    N = np.logspace(np.log10(lo), np.log10(hi), 200)

    curve_rows = component.fabrication_cost_curve(N)
    current_rows = component.fabrication_cost_curve(current_qty)

    fig = _new_figure(figsize=(6.4, 4.2), dpi=150)
    ax = fig.subplots()
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

    return _save_png(fig, path, facecolor=fig.get_facecolor())


def _fmt_input(value) -> str:
    if isinstance(value, str):
        return value
    if value == 0:
        return "0"
    a = abs(value)
    if a >= 1e8 or a < 1e-2:
        return f"{value:.3g}"
    if a >= 1e3:
        return f"{value:,.0f}"
    return f"{value:.4g}"


def tornado_png(result, top_n: int = 15, path: Optional[Union[str, Path]] = None) -> bytes:
    """
    Tornado chart for a sensitivity.SensitivityResult (see
    sensitivity.one_at_a_time).

    One row per input, largest swing at the top, bars drawn from the base
    value to the output at the input's low value (blue) and high value
    (orange). Each row label carries the input's low - high range. Inputs
    with zero swing are left off and counted in the footer, which also
    states the basis (material, processes, volume, quantity, metric, range,
    price source).
    """
    rows = [r for r in result.rows if r["swing"] > 0][:top_n]
    n_zero = sum(1 for r in result.rows if r["swing"] == 0)
    n_hidden = sum(1 for r in result.rows if r["swing"] > 0) - len(rows)
    base = result.base_value
    money = result.units.startswith("$")
    value_fmt = (lambda v: f"${v:,.2f}") if money else (lambda v: f"{v:,.4g}")

    height = 1.6 + 0.38 * max(len(rows), 1)
    fig = _new_figure(figsize=(8.6, height), dpi=150)
    ax = fig.subplots()
    fig.patch.set_facecolor(_SURFACE)
    ax.set_facecolor(_SURFACE)

    low_color, high_color = _CATEGORICAL_PALETTE[0], _CATEGORICAL_PALETTE[1]
    labels = []
    for y, r in enumerate(rows):
        unit = "" if r["units"] in ("-", "level", "level value", "fraction") else f" {r['units']}"
        labels.append(f"{r['label']}  ({_fmt_input(r['low_input'])} – {_fmt_input(r['high_input'])}{unit})")
        for side, color in (("low", low_color), ("high", high_color)):
            out = r[f"{side}_output"]
            if out is None or out == base:
                continue
            ax.barh(y, out - base, left=base, height=0.62, color=color, edgecolor="none", zorder=3,
                    label="Input at low value" if side == "low" else "Input at high value")
            ax.text(out, y, f" {value_fmt(out)} ", ha="left" if out > base else "right", va="center", fontsize=7.5,
                    color=_INK_SECONDARY, zorder=4)

    ax.axvline(base, color=_INK_PRIMARY, linewidth=1, zorder=5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=8, color=_INK_PRIMARY)
    ax.set_ylim(len(rows) - 0.5, -0.5)  # largest swing on top

    # Leave room for the end-of-bar value labels.
    outs = [o for r in rows for o in (r["low_output"], r["high_output"]) if o is not None] + [base]
    lo, hi = min(outs), max(outs)
    pad = (hi - lo) * 0.18 or abs(base) * 0.05 or 1.0
    ax.set_xlim(lo - pad, hi + pad)

    ax.set_xlabel(f"{result.metric} ({result.units})", color=_INK_SECONDARY, fontsize=9)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: value_fmt(v) if not money else f"${v:,.0f}"))
    ax.tick_params(axis="x", colors=_INK_MUTED, labelsize=8)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color=_GRIDLINE, linewidth=1, zorder=0)
    ax.set_axisbelow(True)
    for spine_name in ("top", "right", "left"):
        ax.spines[spine_name].set_visible(False)
    ax.spines["bottom"].set_color(_BASELINE)

    ax.set_title(f"Base {result.metric.lower()}: {value_fmt(base)}", loc="left", color=_INK_PRIMARY, fontsize=10, pad=8)

    handles, legend_labels = ax.get_legend_handles_labels()
    unique = dict(zip(legend_labels, handles))
    order = [k for k in ("Input at low value", "Input at high value") if k in unique]
    if order:
        ax.legend([unique[k] for k in order], order, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2,
                  frameon=False, fontsize=8, labelcolor=_INK_SECONDARY, handlelength=1.2, borderaxespad=0.2)

    b = result.basis
    footer = [
        f"{b['material']} | {b['processes']} | volume {b['volume_mm3']:,.0f} mm³ | N = {b['production_qty']:g} | "
        f"range {b['default_range']}",
        f"Overrides: {b['overrides']}. Element prices: {b['element_prices']}.",
    ]
    extra = []
    if n_hidden:
        extra.append(f"{n_hidden} smaller inputs not shown (top {top_n})")
    if n_zero:
        extra.append(f"{n_zero} inputs with zero swing not shown")
    if extra:
        footer.append("; ".join(extra) + ". Full table: result.rows / to_csv().")

    fig.tight_layout()
    footer_text = fig.text(0.01, 0.0, "\n".join(footer), ha="left", va="top", fontsize=6.5, color=_INK_MUTED,
                           transform=fig.transFigure)

    return _save_png(fig, path, facecolor=fig.get_facecolor(), bbox_extra_artists=(footer_text,), bbox_inches="tight")
