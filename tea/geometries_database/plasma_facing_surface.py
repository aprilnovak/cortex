"""
Plasma Facing Surface


########################################################################
DESCRIPTION
########################################################################
The current pre-defined component geometry. Coefficients (Cc, Cs, Ct,
Cf) are 1.0 (ideal design) for every currently defined process:

    * Cc = 1.0 for shape class 'C1' (basic features, flat/thin wall
      section) - matches the DFM cost reference's Cc charts, which all
      read 1.0 at classification band 1 for every process.
    * Cs = 1.0 at a 2 mm section thickness. Only verifiable against the
      reference for CNC (its machining Cs chart reads ~1.0 by ~2 mm
      section); the reference doesn't cover Hot Rolling, Cold Rolling,
      HIP, Electron Beam, Diffusion Bonding, or Spray Deposition at all,
      so 1.0 is a placeholder default for those until better data exists.
    * Ct = 1.0 at a 0.1 mm tolerance (the value used here, rather than
      the reference's own ideal-design tolerance, so this preset stays
      internally consistent: this geometry's own attributes ARE its
      'ideal' point).
    * Cf = 1.0 at a 1 um Ra surface finish - matches the CNC surface
      finish chart, which reads ~1.0 at Ra = 1 um (the right edge of its
      plotted domain).

########################################################################
PARAMETERS
########################################################################
------------------------------------------------------------------------
name : str
    Name of the geometry. Must be unique

------------------------------------------------------------------------
volume_mm3 : float
    Component volume (mm^3)

------------------------------------------------------------------------
section_thickness_mm : float
    Representative section/wall thickness (mm)

------------------------------------------------------------------------
shape_class : str
    Shape classification band (e.g. 'C1')

------------------------------------------------------------------------
tolerance_mm : float
    Tightest critical total tolerance (mm)

------------------------------------------------------------------------
surface_finish_um_ra : float
    Tightest critical surface finish (um Ra)

------------------------------------------------------------------------
Cc_map, Cs_map, Ct_map, Cf_map : dict
    Dictionaries mapping process name to the relative cost coefficient
    for that process, given this geometry.
"""

name = 'Plasma Facing Surface'

volume_mm3 = 3000000

section_thickness_mm = 2

shape_class = 'C1'

tolerance_mm = 0.1

surface_finish_um_ra = 1

Cc_map = {
    'Hot Rolling': 1.0,
    'Cold Rolling': 1.0,
    'HIP': 1.0,
    'CNC': 1.0,
    'Electron Beam': 1.0,
    'Diffusion Bonding': 1.0,
    'Spray Deposition': 1.0,
}

Cs_map = {
    'Hot Rolling': 1.0,
    'Cold Rolling': 1.0,
    'HIP': 1.0,
    'CNC': 1.0,
    'Electron Beam': 1.0,
    'Diffusion Bonding': 1.0,
    'Spray Deposition': 1.0,
}

Ct_map = {
    'Hot Rolling': 1.0,
    'Cold Rolling': 1.0,
    'HIP': 1.0,
    'CNC': 1.0,
    'Electron Beam': 1.0,
    'Diffusion Bonding': 1.0,
    'Spray Deposition': 1.0,
}

Cf_map = {
    'Hot Rolling': 1.0,
    'Cold Rolling': 1.0,
    'HIP': 1.0,
    'CNC': 1.0,
    'Electron Beam': 1.0,
    'Diffusion Bonding': 1.0,
    'Spray Deposition': 1.0,
}
