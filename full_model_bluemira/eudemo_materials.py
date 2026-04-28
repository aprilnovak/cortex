"""
eudemo_materials.py
===================
Volume-fraction recipes for WCLL, HCLL, and HCPB breeder blankets.
Materials are passed in from cfg (inputs.py) — no construction here.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import List
import sys

import openmc
import inputs as cfg
import materials

# Breeder type enum 
class BreederType(str, Enum):
    WCLL = "WCLL"
    HCLL = "HCLL"
    HCPB = "HCPB"

# Structural / coil / shielding materials 
# Identical regardless of breeder concept.
#ss316    = materials.ss316Ln_ig(7.93)
ccz      = materials.CuCrZr(8.9)
he_low   = materials.Helium(0.0001785)
nb3sn    = materials.Nb3Sn(5.7)
epoxy    = materials.Epoxy(1.207)
bronze   = materials.Bronze(8.8775)
nbti     = materials.NbTi(6.538)
cu       = materials.Cu(8.96)
ss304_b4 = materials.ss304_b4(7.8)

# Plasma fill 
plasma_mat = openmc.Material(name="plasma")
plasma_mat.set_density("g/cm3", 1e-6)
plasma_mat.add_element("H", 1.0)


# Recipe selector
def build_mix_recipes(cfg) -> dict:
    breeder_type = BreederType(cfg.BREEDER_TYPE)

    breeder    = cfg.breeder_material
    coolant    = cfg.coolant_material
    structural = cfg.structural_material
    multiplier = cfg.multiplier_material   # None for WCLL / HCLL
    armor      = cfg.armor_material
    vv         = cfg.vv_material

    shared: dict = {
        "Divertor":            {ccz: 0.00552, cu: 0.00438, structural: 0.5238,
                                armor: 0.01026, coolant: 0.45604},
        "VV_IB":               {vv: 0.6, coolant: 0.4},
        "VV_OB":               {vv: 0.6, coolant: 0.4},
        "VV_ports_all":        {vv: 0.6, coolant: 0.4},
        "VV_port_Fill":        {vv: 0.6, coolant: 0.4},
        "PC_PFC":              {nbti: 0.02895, cu: 0.1169, epoxy: 0.18,
                                bronze: 0.0735, he_low: 0.1682, vv: 0.43245},
        "TFcoil":              {nb3sn: 0.02895, cu: 0.1169, epoxy: 0.18,
                                bronze: 0.0735, he_low: 0.1682, vv: 0.43245},
        "Plasma_Region":       {plasma_mat: 1.00},
        "Cryostat":            {vv: 1.00},
        "RadiationShield_all": {ss304_b4: 1.00},
    }

    layer_builders = {
        BreederType.WCLL: _wcll_layers,
        BreederType.HCLL: _hcll_layers,
        BreederType.HCPB: _hcpb_layers,
    }

    blanket = layer_builders[breeder_type](armor, breeder, coolant, structural, multiplier)

    full = {**blanket, **shared}

    if cfg.SIM_TYPE == "fast_slab":
        # Only keep materials present in the stripped equatorial OB chunk
        _fast_slab_keep = {"Armor", "First_Wall", "VV_OB"}
        # Also keep all OB_Layer_* keys
        full = {
            k: v for k, v in full.items()
            if k in _fast_slab_keep or k.startswith("OB_Layer_")
        }

    return full



def _wcll_layers(armor, breeder, coolant, structural, multiplier) -> dict:
    ib = {
        "Armor":       {armor: 1.00},
        "First_Wall":  {coolant: 0.143, structural: 0.857},
        "IB_Layer_1":  {breeder: 0.833, coolant: 0.028, structural: 0.139},
        "IB_Layer_2":  {breeder: 0.858, coolant: 0.018, structural: 0.124},
        "IB_Layer_3":  {breeder: 0.813, coolant: 0.016, structural: 0.171},
        "IB_Layer_4":  {breeder: 0.813, coolant: 0.016, structural: 0.171},
        "IB_Layer_5":  {breeder: 0.813, coolant: 0.016, structural: 0.171},
        "IB_Layer_6":  {breeder: 0.662, coolant: 0.016, structural: 0.322},
        "IB_Layer_7":  {breeder: 0.586, coolant: 0.016, structural: 0.398},
        "IB_Layer_8":  {coolant: 0.016, structural: 0.984},
        "IB_Layer_9":  {coolant: 0.857, structural: 0.143},
        "IB_Layer_10": {coolant: 0.016, structural: 0.984},
    }
    ob = {k.replace("IB_", "OB_"): v for k, v in ib.items() if k.startswith("IB_")}
    return {**ib, **ob}


def _hcll_layers(armor, breeder, coolant, structural, multiplier) -> dict:
    ib = {
        "Armor":       {armor: 1.00},
        "First_Wall":  {coolant: 0.049, structural: 0.951},
        "IB_Layer_1":  {breeder: 0.883, coolant: 0.002, structural: 0.115},
        "IB_Layer_2":  {breeder: 0.878, coolant: 0.003, structural: 0.119},
        "IB_Layer_3":  {breeder: 0.878, coolant: 0.003, structural: 0.119},
        "IB_Layer_4":  {breeder: 0.313, coolant: 0.118, structural: 0.569},
        "IB_Layer_5":  {breeder: 0.901, coolant: 0.009, structural: 0.090},
        "IB_Layer_6":  {breeder: 0.104, coolant: 0.058, structural: 0.838},
        "IB_Layer_7":  {breeder: 0.010, coolant: 0.861, structural: 0.129},
        "IB_Layer_8":  {breeder: 0.010, coolant: 0.876, structural: 0.114},
        "IB_Layer_9":  {breeder: 0.010, coolant: 0.011, structural: 0.979},
        "IB_Layer_10": {breeder: 0.010, coolant: 0.011, structural: 0.979},
    }
    ob = {k.replace("IB_", "OB_"): v for k, v in ib.items() if k.startswith("IB_")}
    return {**ib, **ob}


def _hcpb_layers(armor, breeder, coolant, structural, multiplier) -> dict:
    ib = {
        "Armor":       {armor: 1.00},
        "First_Wall":  {coolant: 0.282, structural: 0.718},
        "IB_Layer_1":  {breeder: 0.011, multiplier: 0.339, coolant: 0.567, structural: 0.083},
        "IB_Layer_2":  {breeder: 0.054, multiplier: 0.339, coolant: 0.528, structural: 0.079},
        "IB_Layer_3":  {breeder: 0.054, multiplier: 0.339, coolant: 0.530, structural: 0.077},
        "IB_Layer_4":  {breeder: 0.054, multiplier: 0.339, coolant: 0.530, structural: 0.077},
        "IB_Layer_5":  {breeder: 0.054, multiplier: 0.339, coolant: 0.530, structural: 0.077},
        "IB_Layer_6":  {breeder: 0.055, multiplier: 0.028, coolant: 0.643, structural: 0.274},
        "IB_Layer_7":  {coolant: 0.791, structural: 0.209},
        "IB_Layer_8":  {coolant: 0.791, structural: 0.209},
        "IB_Layer_9":  {coolant: 0.040, structural: 0.960},
        "IB_Layer_10": {coolant: 0.040, structural: 0.960},
    }
    ob = {k.replace("IB_", "OB_"): v for k, v in ib.items() if k.startswith("IB_")}
    return {**ib, **ob}