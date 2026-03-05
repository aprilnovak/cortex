# ======================================================
# FUSION 360 PYTHON API SCRIPT 
# TO GENERATE AN EDITED TOKAMAK MODEL
# BASED ON BLUEMIRA OUTPUT WITH:
# 1) IB radial and poloidal discretization 
# 2) OB radial and poloidal discretization
# 3) VV poloidal discretization
# 4) Introduction of Armor and First Wall layers
# ======================================================
# Main steps:
# (0.a) Introduction of user inputs
# (0.b) Function definitions
# (1) Import original Bluemira .step file
# (2) Set units to [cm]
# (3) Create toroidal plasma region for SDR/RandomRay
# (4) Create wedge planes (0 and 22.5 degrees)
# (5) Split bodies by wedge planes.
# (6) Remove bodies outisde of the wedge (keep only between 0-22.5 degrees) (1/16 reactor)
# (7) Delete occurrences (poloidal supporting structures and thermal shield)
# (8) Rename non-breeder components
# (9) Merge Poloidal and Toroidal Coils (reduce complex surfaces for meshing)
# (10) Vacuum Vessel operations
# (10.1) Split ports from VV and rename them
# (10.2) Move ports away from VV (about 2-5 cm -> avoid meshing issues) 
# (10.3) VV face cleanup (remove small faces and heal)
# (10.4) Boolean ports vs VV (will be removed since 10.2 introduced)
# (11) Blanket layer (radial) discretization
# (12) Create IB and OB planes for poloidal discretization
# (13) Split OB, IB, and VV volumes by planes generated in (12)
# (14) Finalize breeders (sort OB and IB chunks with: Armor, First Wall, OB/IB_layers, VV)
# (15) Finalize geometry:
# (15.1) Add Port Filling for SDR
# (15.2) Boolean cuts to avoid overlap (blankets and VV)
# (15.3) Perform overlap checks
# (16) Export STEP file  + Tokamak_inputs.json (for Cubit/OpenMC input)
# (17) Output messages
# ======================================================

import adsk.core
import adsk.fusion
import traceback
import math
import os
import re
import time
import json
from datetime import datetime

# -----------------------------
# USER INPUTS
# -----------------------------
STEP_PATH = r"C:\Users\tpaga\OneDrive\Desktop\UIUC\step_files\EUDEMO.stp"
export_folder = r"C:\Users\tpaga\OneDrive\Desktop\UIUC\step_files\UIUC"

# -----------------------------
# Feel free to change these ones:
# Plasma region
PLASMA_MAJOR_RADIUS_CM = 800.0
PLASMA_TUBE_RADIUS_CM  = 200.0

# First layers offsets
ARMOR_CM = +2.0
FW_CM    = +1.8

#  OB/IB offsets (Still dont have many safe guards against this length)
OB_OFFSETS_CM = [-8.0, -18.0, -28.0, -38.0, -48.0, -58.0, -68.0] # equatorial limit is 95-100 cm (minimal region 80 ish)
IB_OFFSETS_CM = [-8.0, -18.0, -28.0, -38.0, -48.0, -58.0, -68.0] # equatiorial limit around 75-77 cm (minimal 75-77)
#FRACTIONS = [0.125 * i for i in range(1, 8)]  # 0.125..0.875
FRACTIONS = [0.1 * i for i in range(1, 10)]  # 0.125..0.875

RENAME_RULES = [
    ("divertor_mat_homogenised_divertor_2015", "Divertor"),
    ("tfcoil_mat_toroidal_field_coil_2015",    "TFcoil"),
    ("poloidal_coils_mat_poloidal_field_coil", "PFC"),
    ("cryostat_mat_ss316_ln",                  "CR"),
    ("radiationshield_mat_ss316_ln",           "RS1"),
    ("radiationshield",                        "RS2"),
]

# -----------------------------
# Avoid touching these:
# -----------------------------
ANGLE_MAX_DEG = 22.5
ANG_TOL_DEG = 1e-6

REMOVE_NAME_SUBSTRINGS = [
    "coil_structures_mat_ss316_ln",
    "thermal_shield",
]

OB_HINT = "Blanket_mat_Homogenised_HCPB_2015_v3_OB"
IB_HINT = "Blanket_mat_Homogenised_HCPB_2015_v3_IB"
VV_HINT = "vacuumvessel_mat_ss316_ln"

FINAL_COMP_NAME_HINT = "BLANKET_FINAL"
UNITS_SCALE = 1.0
DELETE_ORIGINAL_BLANKET_OCCURRENCES = True
SET_DISPLAY_UNITS_TO_CM = True

OB_NTH_FACE = 4
IB_NTH_FACE = 4

SIDE_CHOICE = "outer"
SPLINE_SAMPLE_POINTS = 140
EDGE_TRIES = 25
ZSPAN_WEIGHT = 2000.0
LEN_WEIGHT   = 1.0
OB_PLANE_NAME_FMT = "{prefix}_PLANE_{frac:04d}"
IB_PLANE_NAME_FMT = "{prefix}_PLANE_{frac:04d}"

FINALIZE_DELETE_ORIGINALS = False
FINALIZE_HIDE_ORIGINALS   = True

VV_KEEP_FACES    = 4
VV_MIN_FACE_AREA = 0.0
VV_HEAL_ON_DELETE = True

# -----------------------------
# Toggles
# -----------------------------
DO_CREATE_PLASMA_REGION = False  

# -----------------------------
# HEAVY STEP TOGGLES
# -----------------------------
DO_PFC_TRIPLET_MERGE   = True
DO_TFCOIL_JOIN         = True

DO_VV_OPS_EARLY        = True
DO_PORT_FILLING        = True

DO_DEOVERLAP_15_2       = True
DO_BF_X_VV_CUT_15_3    = True
DO_CUT_DIVERTORS       = True
DO_OVERLAP_VERIFY      = True

# Port filling correction:
PORT_FILL_SWAP_WH = True # function swapped H/L (fix tbd later)
PORT_FILL_PORTNAME = "VV_port_2"
PORT_FILL_THETA_DEG = 11.25
PORT_FILL_WIDTH_CM  = 180.0
PORT_FILL_HEIGHT_CM = 250.0
PORT_FILL_USE_SQUARE = False
PORT_FILL_DEPTH_SCALE = 0.60

# ======================================================
# Small utilities
# ======================================================
def _safe_name(x) -> str:
    try:
        return (x.name or "")
    except Exception:
        return ""

def _try_set_name(x, nm: str):
    try:
        x.name = nm
    except Exception:
        pass

def _log_to_text_palette(app, msg: str):
    try:
        pal = app.userInterface.palettes.itemById('TextCommands')
        if pal:
            pal.writeText(msg)
    except Exception:
        pass

def _collect_occ_tree_from_root(root_comp: adsk.fusion.Component):
    out, stack = [], []
    for i in range(root_comp.occurrences.count):
        stack.append(root_comp.occurrences.item(i))
    while stack:
        occ = stack.pop()
        out.append(occ)
        try:
            for j in range(occ.childOccurrences.count):
                stack.append(occ.childOccurrences.item(j))
        except Exception:
            pass
    return out

def _occ_blob(occ: adsk.fusion.Occurrence) -> str:
    n1 = _safe_name(occ)
    n2 = ""
    try:
        n2 = _safe_name(occ.component)
    except Exception:
        pass
    return (n1 + " " + n2).lower()

def _find_all_occurrences_by_hint(root_comp: adsk.fusion.Component, hint: str):
    h = (hint or "").lower()
    matches = []
    for occ in _collect_occ_tree_from_root(root_comp):
        if h in _occ_blob(occ):
            matches.append(occ)
    return matches

def _find_first_occurrence_by_hint(root_comp: adsk.fusion.Component, hint: str):
    ms = _find_all_occurrences_by_hint(root_comp, hint)
    return ms[0] if ms else None

def _name_matches_remove_list(occ: adsk.fusion.Occurrence) -> bool:
    blob = _occ_blob(occ)
    for s in REMOVE_NAME_SUBSTRINGS:
        if (s or "").lower() in blob:
            return True
    return False

def _solid_proxy_bodies(occ: adsk.fusion.Occurrence):
    try:
        for b in list(occ.bRepBodies):
            try:
                if b.isSolid:
                    yield b
            except Exception:
                continue
    except Exception:
        return

def _solid_definition_bodies(occ: adsk.fusion.Occurrence):
    out = []
    try:
        comp = occ.component
        for i in range(comp.bRepBodies.count):
            b = comp.bRepBodies.item(i)
            try:
                if b.isSolid:
                    out.append(b)
            except Exception:
                pass
    except Exception:
        pass
    return out

def _native_token(obj):
    try:
        o = obj.nativeObject if getattr(obj, "nativeObject", None) else obj
    except Exception:
        o = obj
    try:
        return o.entityToken
    except Exception:
        return None

def _proxy_volume(b: adsk.fusion.BRepBody) -> float:
    try:
        return float(b.physicalProperties.volume)
    except Exception:
        return 0.0

def _rename_body_native_if_possible(body_proxy: adsk.fusion.BRepBody, new_name: str):
    try:
        nb = body_proxy.nativeObject if getattr(body_proxy, "nativeObject", None) else body_proxy
    except Exception:
        nb = body_proxy
    try:
        nb.name = new_name
        return
    except Exception:
        pass
    try:
        body_proxy.name = new_name
    except Exception:
        pass

def _find_body_by_exact_name_in_occ(occ: adsk.fusion.Occurrence, want: str):
    w = (want or "").strip().lower()
    for b in _solid_proxy_bodies(occ):
        if _safe_name(b).strip().lower() == w:
            return b
    return None

def _find_bodies_by_prefix_in_occ(occ: adsk.fusion.Occurrence, prefix: str):
    p = (prefix or "").strip().lower()
    out = []
    for b in _solid_proxy_bodies(occ):
        nm = _safe_name(b).strip().lower()
        if nm.startswith(p):
            out.append(b)
    return out

def _as_native(obj):
    cur = obj
    try:
        while hasattr(cur, "nativeObject") and cur.nativeObject:
            nxt = cur.nativeObject
            if nxt is cur:
                break
            cur = nxt
    except Exception:
        pass
    return cur


class StepTimer:
    def __init__(self):
        self.t0 = time.time()
        self.last = self.t0
        self.rows = []

    def mark(self, label: str):
        now = time.time()
        dt = now - self.last
        self.last = now
        self.rows.append((label, dt))
        return dt

    def total(self):
        return time.time() - self.t0

    def format_table(self) -> str:
        if not self.rows:
            return "(no step timings)"
        w = max(len(lbl) for (lbl, _) in self.rows)
        lines = []
        for lbl, sec in self.rows:
            lines.append(f"{lbl:<{w}} : {sec:8.2f} s")
        lines.append("-" * (w + 14))
        lines.append(f"{'TOTAL':<{w}} : {self.total():8.2f} s")
        return "\n".join(lines)

def _ensure_parent_folder_exists(file_path: str) -> str:
    """
    Ensures parent folder exists for a given file path.
    Returns normalized absolute path.
    """
    try:
        file_path = os.path.abspath(file_path)
        parent = os.path.dirname(file_path)
        if parent and (not os.path.isdir(parent)):
            os.makedirs(parent, exist_ok=True)
        return file_path
    except Exception:
        # if anything weird happens, just return original
        return file_path

# ======================================================
# Split helpers
# ======================================================
def _split_body_by_tool(root_comp: adsk.fusion.Component,
                        target_body_proxy: adsk.fusion.BRepBody,
                        splitting_tool,
                        extend=True):
    split_feats = root_comp.features.splitBodyFeatures
    inp = split_feats.createInput(target_body_proxy, splitting_tool, extend)
    return split_feats.add(inp)

def _plane_frac_from_name(p):
    nm = _safe_name(p)
    m = re.search(r"_PLANE_(\d{4})$", nm.upper())
    return int(m.group(1)) if m else -1

def _plane_in_body_context(pl_root, body_proxy):
    """
    Return pl_root as a proxy plane in the SAME assembly context as body_proxy.
    """
    try:
        occ = body_proxy.assemblyContext
    except Exception:
        occ = None

    pln = _as_native(pl_root)
    if occ:
        try:
            return pln.createForAssemblyContext(occ)
        except Exception:
            return pln
    return pln

def _vec_from_points(p, q):
    return adsk.core.Vector3D.create(q.x - p.x, q.y - p.y, q.z - p.z)

def _try_face_point_and_normal(face: adsk.fusion.BRepFace):
    """
    Returns (Point3D, Vector3D) or (None,None) if cannot compute.
    Uses face.pointOnFace and SurfaceEvaluator.getParametersAtPoints + getNormalAtParameter.
    """
    try:
        p = face.pointOnFace
    except Exception:
        return None, None

    try:
        ev = face.evaluator
    except Exception:
        return p, None

    try:
        pts = adsk.core.ObjectCollection.create()
        pts.add(p)
        ok, params = ev.getParametersAtPoints(pts)
        if not ok or not params or params.count < 1:
            return p, None
        uv = params.item(0)
        u = uv.x
        v = uv.y
    except Exception:
        return p, None

    try:
        ok2, n = ev.getNormalAtParameter(u, v)
        if not ok2:
            return p, None
        return p, n
    except Exception:
        return p, None

def _sort_faces_by_area(face_list, reverse=True):
    scored = []
    for f in face_list:
        scored.append((_face_area(f), f))
    scored.sort(key=lambda t: t[0], reverse=reverse)
    return scored

# ======================================================
# VV operations
# ======================================================
def _largest_face_any(body: adsk.fusion.BRepBody):
    best_f = None
    best_a = -1.0
    for i in range(body.faces.count):
        f = body.faces.item(i)
        try:
            a = float(f.area)
        except Exception:
            continue
        if a > best_a:
            best_a = a
            best_f = f
    return best_f, best_a

def split_and_name_vv_ports(root_comp: adsk.fusion.Component,
                            vv_hint: str = "vacuumvessel_mat_ss316_ln"):
    vv_occ = _find_first_occurrence_by_hint(root_comp, vv_hint)
    if vv_occ is None:
        return (False, 0, "VV: occurrence not found.")

    vv_bodies = list(_solid_proxy_bodies(vv_occ))
    if not vv_bodies:
        return (False, 0, "VV: no solid bodies found.")

    if len(vv_bodies) > 1:
        solids_sorted = sorted(vv_bodies, key=_proxy_volume, reverse=True)
        _rename_body_native_if_possible(solids_sorted[0], "VV_1")
        for i, b in enumerate(solids_sorted[1:], start=1):
            _rename_body_native_if_possible(b, f"VV_port_{i}")
        return (False, len(solids_sorted) - 1, f"VV: already {len(solids_sorted)} bodies, renamed ports.")

    vv_body = vv_bodies[0]
    big_face, big_area = _largest_face_any(vv_body)
    if big_face is None:
        return (False, 0, "VV: could not find largest face.")

    try:
        feat = _split_body_by_tool(root_comp, vv_body, big_face, extend=True)
    except Exception as e:
        return (False, 0, f"VV: split failed: {e}")

    new_bodies = []
    try:
        for i in range(feat.bodies.count):
            nb = feat.bodies.item(i)
            try:
                if nb.isSolid:
                    new_bodies.append(nb)
            except Exception:
                pass
    except Exception:
        pass

    if not new_bodies:
        new_bodies = list(_solid_proxy_bodies(vv_occ))

    if len(new_bodies) < 2:
        _rename_body_native_if_possible(new_bodies[0], "VV_1")
        return (False, 0, f"VV: split area={big_area:.6g} produced 1 body (no ports).")

    solids_sorted = sorted(new_bodies, key=_proxy_volume, reverse=True)
    _rename_body_native_if_possible(solids_sorted[0], "VV_1")
    for i, b in enumerate(solids_sorted[1:], start=1):
        _rename_body_native_if_possible(b, f"VV_port_{i}")

    return (True, len(solids_sorted) - 1,
            f"VV: split area={big_area:.6g}, bodies={len(solids_sorted)}. Named VV_1 + {len(solids_sorted)-1} ports.")

def _face_area(face: adsk.fusion.BRepFace) -> float:
    try:
        return float(face.area)
    except Exception:
        return 0.0

def _delete_faces_feature(comp: adsk.fusion.Component,
                          faces: adsk.core.ObjectCollection,
                          heal: bool = True):
    dff = comp.features.deleteFaceFeatures
    if hasattr(dff, "createInput"):
        inp = dff.createInput(faces, heal)
        return dff.add(inp)
    try:
        return dff.add(faces, heal)
    except TypeError:
        pass
    return dff.add(faces)

def vv1_face_cleanup_keep_largest(vv_occ: adsk.fusion.Occurrence,
                                  keep_faces: int = 4,
                                  min_face_area: float = 0.0,
                                  heal: bool = True):
    solids = _solid_definition_bodies(vv_occ)
    if not solids:
        return (False, "VV1 face-cleanup: no solid definition bodies found.")

    vv1_def = None
    for b in solids:
        if _safe_name(b).strip().lower() == "vv_1":
            vv1_def = b
            break
    if vv1_def is None:
        vv1_def = sorted(solids, key=lambda bb: float(getattr(bb.physicalProperties, "volume", 0.0)), reverse=True)[0]

    face_list = []
    for i in range(vv1_def.faces.count):
        f = vv1_def.faces.item(i)
        face_list.append((_face_area(f), f))
    face_list.sort(key=lambda t: t[0], reverse=True)

    total = len(face_list)
    if total <= keep_faces:
        return (True, f"VV1 face-cleanup: VV_1 has {total} faces (<= keep {keep_faces}), nothing deleted.")

    to_delete = face_list[keep_faces:]
    if min_face_area and min_face_area > 0.0:
        to_delete = [(a, f) for (a, f) in to_delete if a >= min_face_area]
    if not to_delete:
        return (True, f"VV1 face-cleanup: nothing matched delete filter; faces={total}, kept={keep_faces}.")

    fc = adsk.core.ObjectCollection.create()
    for _, f in to_delete:
        fc.add(f)

    _delete_faces_feature(vv_occ.component, fc, heal=heal)
    return (True, f"VV1 face-cleanup: deleted {len(to_delete)} face(s), kept top {keep_faces} (total was {total}).")

def vv_ports_cut_against_vv1(root_comp: adsk.fusion.Component,
                             vv_occ: adsk.fusion.Occurrence,
                             vv1_name: str = "VV_1",
                             port_prefix: str = "VV_port_",
                             keep_tool_vv1: bool = True):
    vv1 = _find_body_by_exact_name_in_occ(vv_occ, vv1_name)
    if vv1 is None:
        return (False, f"Ports vs VV1 cut: '{vv1_name}' not found.")

    ports = []
    for b in list(_solid_proxy_bodies(vv_occ)):
        nm = _safe_name(b).strip().lower()
        if nm.startswith(port_prefix.lower()):
            ports.append(b)

    if not ports:
        return (True, "Ports vs VV1 cut: no ports found.")

    cmb = root_comp.features.combineFeatures
    tools = adsk.core.ObjectCollection.create()
    tools.add(vv1)

    ok = 0
    fail = 0
    first_err = None

    for p in ports:
        try:
            ci = cmb.createInput(p, tools)
            ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
            ci.isKeepToolBodies = bool(keep_tool_vv1)
            cmb.add(ci)
            ok += 1
        except Exception as e:
            fail += 1
            if first_err is None:
                first_err = str(e)

    if fail == 0:
        return (True, f"Ports vs VV1 cut: CUT applied to {ok} port(s) using VV_1 tool (VV_1 unchanged).")
    return (False, f"Ports vs VV1 cut: ok={ok}, fail={fail}. First error: {first_err}")

def split_vv1_by_ob_ib_planes(root_comp: adsk.fusion.Component,
                             vv_occ: adsk.fusion.Occurrence,
                             vv1_name: str,
                             path_data_by_j,
                             path_data_ib_by_j,
                             split_feats_root):
    """
    Split VV_1 using the OB and IB plane families (poloidal discretization),
    in the *correct assembly context*.

    After each split, we re-pick the largest remaining solid and rename it VV_1
    so subsequent planes keep splitting the main vessel rather than tiny scraps.
    """
    vv1 = _find_body_by_exact_name_in_occ(vv_occ, vv1_name)
    if vv1 is None:
        # fallback: pick largest solid in vv_occ
        solids = list(_solid_proxy_bodies(vv_occ))
        if not solids:
            return (False, "VV_1 plane split: no VV solids found.")
        vv1 = sorted(solids, key=_proxy_volume, reverse=True)[0]
        _rename_body_native_if_possible(vv1, vv1_name)

    # collect ALL planes (OB + IB) once
    all_planes_root_native = []
    for j in sorted(path_data_by_j.keys()):
        _, _, _, pls = path_data_by_j[j]
        all_planes_root_native.extend([_as_native(p) for p in pls])

    for j in sorted(path_data_ib_by_j.keys()):
        _, _, _, pls = path_data_ib_by_j[j]
        all_planes_root_native.extend([_as_native(p) for p in pls])

    # sort by fraction key so we split "top-down" consistently
    all_planes_root_native = sorted(all_planes_root_native, key=_plane_frac_from_name, reverse=True)

    ok = 0
    skip = 0
    first_err = None

    for pl_root_native in all_planes_root_native:
        # VV_1 reference changes after every split; refresh it each plane
        vv1 = _find_body_by_exact_name_in_occ(vv_occ, vv1_name)
        if vv1 is None:
            solids = list(_solid_proxy_bodies(vv_occ))
            if not solids:
                break
            vv1 = sorted(solids, key=_proxy_volume, reverse=True)[0]
            _rename_body_native_if_possible(vv1, vv1_name)

        try:
            pl_for_vv = _plane_in_body_context(pl_root_native, vv1)
            inp = split_feats_root.createInput(vv1, pl_for_vv, True)
            feat = split_feats_root.add(inp)

            # if Fusion returns new bodies, keep the largest as VV_1 and leave others as pieces
            new_solids = []
            try:
                for ii in range(feat.bodies.count):
                    nb = feat.bodies.item(ii)
                    if getattr(nb, "isSolid", True):
                        new_solids.append(nb)
            except Exception:
                pass

            if len(new_solids) < 2:
                skip += 1
                continue

            largest = sorted(new_solids, key=_proxy_volume, reverse=True)[0]
            _rename_body_native_if_possible(largest, vv1_name)
            ok += 1

        except Exception as e:
            skip += 1
            if first_err is None:
                first_err = str(e)

    msg = f"VV_1 plane split: ok={ok}, skip={skip}"
    if first_err:
        msg += f". First error: {first_err}"
    return (ok > 0), msg

def vv_poloidal_split_directed(root_comp: adsk.fusion.Component,
                              vv_occ: adsk.fusion.Occurrence,
                              ob_occ: adsk.fusion.Occurrence,
                              ib_occ: adsk.fusion.Occurrence):
    try:
        vv1 = _find_body_by_exact_name_in_occ(vv_occ, "VV_1")
        if vv1 is None:
            return (False, "VV poloidal split: VV_1 not found.")

        ob_defs = _solid_definition_bodies(ob_occ)
        ib_defs = _solid_definition_bodies(ib_occ)
        if not ob_defs or not ib_defs:
            return (False, "VV poloidal split: OB/IB breeder bodies not found.")

        ob_main = sorted(ob_defs, key=lambda b: float(getattr(b.physicalProperties, "volume", 0.0)), reverse=True)[0]
        ib_main = sorted(ib_defs, key=lambda b: float(getattr(b.physicalProperties, "volume", 0.0)), reverse=True)[0]

        ob_com = ob_main.physicalProperties.centerOfMass
        ib_com = ib_main.physicalProperties.centerOfMass

        ib_faces = [ib_main.faces.item(i) for i in range(ib_main.faces.count)]
        ob_faces = [ob_main.faces.item(i) for i in range(ob_main.faces.count)]
        if len(ib_faces) < 2 or len(ob_faces) < 1:
            return (False, "VV poloidal split: not enough breeder faces.")

        ib_small_sorted = _sort_faces_by_area(ib_faces, reverse=False)
        ob_small_sorted = _sort_faces_by_area(ob_faces, reverse=False)

        f1 = ib_small_sorted[1][1]  # 2nd smallest IB1 face
        f2 = ib_small_sorted[0][1]  # smallest IB1 face
        f3 = ob_small_sorted[0][1]  # smallest OB2 face (use OB main pre-layer)

        try:
            _split_body_by_tool(root_comp, vv1, f1, extend=True)
        except Exception as e:
            return (False, f"VV poloidal split: first cut failed: {e}")

        bodies_after1 = list(_solid_proxy_bodies(vv_occ))
        candidates = [b for b in bodies_after1 if _safe_name(b).lower().startswith("vv_") and "port" not in _safe_name(b).lower()]
        candidates = sorted(candidates, key=_proxy_volume, reverse=True)
        if len(candidates) < 2:
            candidates = sorted(list(_solid_proxy_bodies(vv_occ)), key=_proxy_volume, reverse=True)
        if len(candidates) < 2:
            return (False, "VV poloidal split: first cut produced <2 solids.")

        a, b = candidates[0], candidates[1]
        a_com = a.physicalProperties.centerOfMass
        b_com = b.physicalProperties.centerOfMass

        def _d2(p, q):
            dx = p.x - q.x
            dy = p.y - q.y
            dz = p.z - q.z
            return dx*dx + dy*dy + dz*dz

        if _d2(a_com, ib_com) <= _d2(b_com, ib_com):
            vv_ib = a
            vv_ob = b
        else:
            vv_ib = b
            vv_ob = a

        _rename_body_native_if_possible(vv_ib, "VV_IB")
        _rename_body_native_if_possible(vv_ob, "VV_OB")

        vv_ib = _find_body_by_exact_name_in_occ(vv_occ, "VV_IB")
        vv_ob = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB")
        if vv_ib is None or vv_ob is None:
            return (False, "VV poloidal split: could not resolve VV_IB/VV_OB after first cut.")

        # Cut 2: VV_IB by IB smallest -> divertor_1
        try:
            _split_body_by_tool(root_comp, vv_ib, f2, extend=True)
        except Exception:
            pass

        vv_div1 = None
        vv_ib_main = None
        ib_candidates = _find_bodies_by_prefix_in_occ(vv_occ, "vv_ib")
        if ib_candidates:
            best = None
            bestd = 1e99
            for bb in ib_candidates:
                cc = bb.physicalProperties.centerOfMass
                d = _d2(cc, ib_com)
                if d < bestd:
                    bestd = d
                    best = bb
            vv_ib_main = best
            rem = [x for x in ib_candidates if x != vv_ib_main]
            if rem:
                vv_div1 = sorted(rem, key=_proxy_volume, reverse=True)[0]

        if vv_ib_main:
            _rename_body_native_if_possible(vv_ib_main, "VV_IB")
        if vv_div1:
            _rename_body_native_if_possible(vv_div1, "VV_divertor_1")

        # Cut 3: VV_OB by OB smallest -> divertor_2
        try:
            _split_body_by_tool(root_comp, vv_ob, f3, extend=True)
        except Exception:
            pass

        vv_div2 = None
        vv_ob_main = None
        ob_candidates = _find_bodies_by_prefix_in_occ(vv_occ, "vv_ob")
        if ob_candidates:
            best = None
            bestd = 1e99
            for bb in ob_candidates:
                cc = bb.physicalProperties.centerOfMass
                d = _d2(cc, ob_com)
                if d < bestd:
                    bestd = d
                    best = bb
            vv_ob_main = best
            rem = [x for x in ob_candidates if x != vv_ob_main]
            if rem:
                vv_div2 = sorted(rem, key=_proxy_volume, reverse=True)[0]

        if vv_ob_main:
            _rename_body_native_if_possible(vv_ob_main, "VV_OB")
        if vv_div2:
            _rename_body_native_if_possible(vv_div2, "VV_divertor_2")

        div1 = _find_body_by_exact_name_in_occ(vv_occ, "VV_divertor_1")
        div2 = _find_body_by_exact_name_in_occ(vv_occ, "VV_divertor_2")

        if div1 and div2:
            cmb = root_comp.features.combineFeatures
            tools = adsk.core.ObjectCollection.create()
            tools.add(div2)
            ci = cmb.createInput(div1, tools)
            ci.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
            ci.isKeepToolBodies = False
            cmb.add(ci)
            _rename_body_native_if_possible(div1, "VV_divertor")
            try:
                _rename_body_native_if_possible(div2, "VV_divertor_old")
            except Exception:
                pass
        elif div1:
            _rename_body_native_if_possible(div1, "VV_divertor")
        elif div2:
            _rename_body_native_if_possible(div2, "VV_divertor")

        solids = list(_solid_proxy_bodies(vv_occ))
        return (True, f"VV poloidal split: created VV_OB, VV_IB + joined divertor parts into VV_divertor. solids={len(solids)}")

    except Exception as e:
        return (False, f"VV poloidal split: exception: {e}")

def _pick_face_by_rank(body: adsk.fusion.BRepBody, rank: int, largest=True):
    faces = [body.faces.item(i) for i in range(body.faces.count)]
    if not faces:
        return None
    scored = _sort_faces_by_area(faces, reverse=largest)
    idx = max(0, min(rank - 1, len(scored) - 1))
    return scored[idx][1]


def vv_split_ob_into_3(root_comp, vv_occ, ob_occ_for_faces):
    vv_ob = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB")
    if vv_ob is None:
        return (False, "VV breeder split (OB): VV_OB not found.")

    ob_defs = _solid_definition_bodies(ob_occ_for_faces)
    if not ob_defs:
        return (False, "VV breeder split (OB): OB breeder bodies not found.")
    ob_main = sorted(ob_defs, key=lambda b: float(getattr(b.physicalProperties, "volume", 0.0)), reverse=True)[0]

    f2 = _pick_face_by_rank(ob_main, 2, largest=True)  # 2nd largest
    f3 = _pick_face_by_rank(ob_main, 3, largest=True)  # 3rd largest
    if f2 is None or f3 is None:
        return (False, "VV breeder split (OB): could not pick 2nd/3rd largest faces.")

    try:
        _split_body_by_tool(root_comp, vv_ob, f2, extend=True)
    except Exception as e:
        return (False, f"VV breeder split (OB): first split failed: {e}")

    ob_pieces = _find_bodies_by_prefix_in_occ(vv_occ, "vv_ob")
    vv_ob = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB")
    if vv_ob is None:
        if not ob_pieces:
            return (False, "VV breeder split (OB): no VV_OB* bodies after first split.")
        vv_ob = sorted(ob_pieces, key=_proxy_volume, reverse=True)[0]
        _rename_body_native_if_possible(vv_ob, "VV_OB")

    try:
        _split_body_by_tool(root_comp, vv_ob, f3, extend=True)
    except Exception as e:
        return (False, f"VV breeder split (OB): second split failed: {e}")

    b0 = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB")
    b1 = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB (1)") or _find_body_by_exact_name_in_occ(vv_occ, "VV_OB(1)")
    b2 = _find_body_by_exact_name_in_occ(vv_occ, "VV_OB (2)") or _find_body_by_exact_name_in_occ(vv_occ, "VV_OB(2)")

    all_ob = _find_bodies_by_prefix_in_occ(vv_occ, "vv_ob")
    uniq = []
    for bb in all_ob:
        if bb not in uniq:
            uniq.append(bb)
    if len(uniq) < 3:
        return (False, f"VV OB rename: found {len(uniq)} VV_OB* bodies, expected 3.")

    if b0 and b1 and b2:
        _rename_body_native_if_possible(b1, "VV_OB1")
        _rename_body_native_if_possible(b0, "VV_OB2")
        _rename_body_native_if_possible(b2, "VV_OB3")
        return (True, "VV breeder split (OB): split into 3 and renamed VV_OB1,VV_OB2,VV_OB3.")
    else:
        base = None
        for bb in uniq:
            if _safe_name(bb).strip().lower() == "vv_ob":
                base = bb
                break
        if base is None:
            base = sorted(uniq, key=_proxy_volume, reverse=True)[0]
        rem = [bb for bb in uniq if bb != base]
        rem_sorted = sorted(rem, key=_proxy_volume, reverse=False)
        _rename_body_native_if_possible(rem_sorted[0], "VV_OB1")
        _rename_body_native_if_possible(base, "VV_OB2")
        _rename_body_native_if_possible(rem_sorted[1], "VV_OB3")
        return (True, "VV breeder split (OB): split into 3 and renamed (fallback mapping).")

def vv_split_ib_into_2_by_facing_face(root_comp, vv_occ, ib_occ1_for_face, ib_occ2_for_direction):
    vv_ib = _find_body_by_exact_name_in_occ(vv_occ, "VV_IB")
    if vv_ib is None:
        return (False, "VV breeder split (IB): 'VV_IB' not found in VV occurrence.")

    ib1_defs = _solid_definition_bodies(ib_occ1_for_face)
    ib2_defs = _solid_definition_bodies(ib_occ2_for_direction)
    if not ib1_defs or not ib2_defs:
        return (False, "VV breeder split (IB): IB breeder bodies not found.")

    ib1 = sorted(ib1_defs, key=lambda b: float(getattr(b.physicalProperties, "volume", 0.0)), reverse=True)[0]
    ib2 = sorted(ib2_defs, key=lambda b: float(getattr(b.physicalProperties, "volume", 0.0)), reverse=True)[0]

    c2 = ib2.physicalProperties.centerOfMass

    best = None
    best_score = -1e99

    for i in range(ib1.faces.count):
        f = ib1.faces.item(i)
        a = _face_area(f)
        if a <= 0.0:
            continue
        p, n = _try_face_point_and_normal(f)
        if p is None or n is None:
            continue
        v = _vec_from_points(p, c2)
        try:
            v.normalize()
            nn = adsk.core.Vector3D.create(n.x, n.y, n.z)
            nn.normalize()
            score = nn.dotProduct(v) * a
        except Exception:
            score = -1e99
        if score > best_score:
            best_score = score
            best = f

    if best is None:
        best = _pick_face_by_rank(ib1, 1, largest=True)

    if best is None:
        return (False, "VV breeder split (IB): could not select a facing face on IB1.")

    try:
        _split_body_by_tool(root_comp, vv_ib, best, extend=True)
    except Exception as e:
        return (False, f"VV breeder split (IB): split failed: {e}")

    b0 = _find_body_by_exact_name_in_occ(vv_occ, "VV_IB")
    b1 = _find_body_by_exact_name_in_occ(vv_occ, "VV_IB (1)") or _find_body_by_exact_name_in_occ(vv_occ, "VV_IB(1)")

    ib_pieces = _find_bodies_by_prefix_in_occ(vv_occ, "vv_ib")
    uniq = []
    for bb in ib_pieces:
        if bb not in uniq:
            uniq.append(bb)
    if len(uniq) < 2:
        return (False, f"VV breeder split (IB): expected 2 VV_IB* bodies, got {len(uniq)}.")

    if b0 and b1:
        _rename_body_native_if_possible(b0, "VV_IB1")
        _rename_body_native_if_possible(b1, "VV_IB2")
        return (True, "VV breeder split (IB): split into 2 and renamed VV_IB1,VV_IB2.")
    else:
        c1 = ib1.physicalProperties.centerOfMass
        def _d2(p,q):
            dx=p.x-q.x; dy=p.y-q.y; dz=p.z-q.z
            return dx*dx+dy*dy+dz*dz
        bestb = None
        bestd = 1e99
        for bb in uniq:
            cc = bb.physicalProperties.centerOfMass
            d = _d2(cc, c1)
            if d < bestd:
                bestd = d
                bestb = bb
        other = [bb for bb in uniq if bb != bestb][0]
        _rename_body_native_if_possible(bestb, "VV_IB1")
        _rename_body_native_if_possible(other, "VV_IB2")
        return (True, "VV breeder split (IB): split into 2 and renamed (fallback mapping).")

# ======================================================
# Component helpers + layering
# ======================================================
def new_component(parent_comp, name: str):
    occ = parent_comp.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    _try_set_name(occ.component, name)
    return occ, occ.component

def copy_body_into_def(comp, src_body, new_name: str):
    feat = comp.features.copyPasteBodies.add(src_body)
    if not feat or feat.bodies.count < 1:
        raise RuntimeError(f"Copy failed into '{_safe_name(comp)}'.")
    b = feat.bodies.item(0)
    _try_set_name(b, new_name)
    return b

def nth_largest_face(body, n: int):
    faces = [body.faces.item(i) for i in range(body.faces.count)]
    if not faces:
        return None
    scored = []
    for f in faces:
        try:
            a = float(f.area)
        except Exception:
            a = 0.0
        scored.append((a, f))
    scored.sort(key=lambda t: t[0], reverse=True)
    idx = max(0, min(n - 1, len(scored) - 1))
    return scored[idx][1]

def offset_body_by_nth_face(comp, body, nth_face: int, offset_val: float):
    seed = nth_largest_face(body, nth_face)
    if seed is None:
        raise RuntimeError(f"No faces found on body '{_safe_name(body)}'.")
    off_feats = comp.features.offsetFacesFeatures
    inp = off_feats.createInput([seed], adsk.core.ValueInput.createByReal(offset_val))
    inp.isChainSelection = True
    off_feats.add(inp)

def cut_shell(comp, outer_tool, inner_tool, result_name: str):
    shell = copy_body_into_def(comp, outer_tool, result_name)
    tools = adsk.core.ObjectCollection.create()
    tools.add(inner_tool)
    cmb = comp.features.combineFeatures
    ci = cmb.createInput(shell, tools)
    ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
    ci.isKeepToolBodies = True
    cmb.add(ci)
    return shell

def build_layers_for_breeder(root_comp: adsk.fusion.Component,
                            parent_final_comp: adsk.fusion.Component,
                            src_bodies,
                            tag: str,
                            nth_face: int,
                            layer_offsets_cm):
    tag_occ, tag_comp = new_component(parent_final_comp, tag)

    neg_cm = sorted([v for v in layer_offsets_cm if v < 0], reverse=True)
    per_body_summary = []

    for bi, src_body in enumerate(src_bodies, start=1):
        build_occ, build_comp = new_component(root_comp, f"{tag}_BUILD_B{bi}")

        base = copy_body_into_def(build_comp, src_body, "Base")

        armor_tool = copy_body_into_def(build_comp, base, "P_Armor")
        offset_body_by_nth_face(build_comp, armor_tool, nth_face, ARMOR_CM * UNITS_SCALE)

        fw_tool = copy_body_into_def(build_comp, base, "P_FW")
        offset_body_by_nth_face(build_comp, fw_tool, nth_face, FW_CM * UNITS_SCALE)

        neg_tools = {}
        for v in neg_cm:
            nm = f"M_{int(abs(v))}"
            tb = copy_body_into_def(build_comp, base, nm)
            offset_body_by_nth_face(build_comp, tb, nth_face, float(v) * UNITS_SCALE)
            neg_tools[v] = tb

        armor_shell = cut_shell(build_comp, armor_tool, fw_tool, "Armor_shell")
        fw_shell    = cut_shell(build_comp, fw_tool, base, "FirstWall_shell")

        layer_bodies = []
        if neg_cm:
            layer_bodies.append(cut_shell(build_comp, base, neg_tools[neg_cm[0]], "Breeder_1"))
            for li in range(1, len(neg_cm)):
                outer = neg_tools[neg_cm[li - 1]]
                inner = neg_tools[neg_cm[li]]
                layer_bodies.append(cut_shell(build_comp, outer, inner, f"Breeder_{li + 1}"))
            core = copy_body_into_def(build_comp, neg_tools[neg_cm[-1]], f"Breeder_{len(neg_cm) + 1}")
            layer_bodies.append(core)
        else:
            layer_bodies.append(copy_body_into_def(build_comp, base, "Breeder_1"))

        bfinal_occ, bfinal_comp = new_component(tag_comp, f"{tag}{bi}")
        copy_body_into_def(bfinal_comp, armor_shell, "Armor_shell")
        copy_body_into_def(bfinal_comp, fw_shell,    "FirstWall_shell")
        for li, body in enumerate(layer_bodies, start=1):
            copy_body_into_def(bfinal_comp, body, f"{tag}_Layer_{li}")

        try:
            build_occ.deleteMe()
        except Exception:
            pass

        per_body_summary.append(f"{tag}{bi}: layers={len(layer_bodies)} (includes core), nth_face={nth_face}")

    return per_body_summary, tag_comp

# ======================================================
# Renaming, Inventory, and JSON
# ======================================================
def rename_blanket_final(blanket_final_comp: adsk.fusion.Component):
    def rename_bodies_in_comp(c, prefix: str):
        for i in range(c.bRepBodies.count):
            bb = c.bRepBodies.item(i)
            old = _safe_name(bb)
            low = old.lower()
            if low == "armor_shell":
                new = f"{prefix}_Armor"
            elif low in ("firstwall_shell", "first_wall_shell", "first wall_shell"):
                new = f"{prefix}_First_Wall"
            elif "layer" in low:
                nums = re.findall(r"\d+", old)
                new = f"{prefix}_Layer_{nums[-1]}" if nums else f"{prefix}_Layer"
            else:
                safe = old.replace(" ", "_")
                new = f"{prefix}_{safe}" if safe else f"{prefix}_Body_{i+1}"
            _try_set_name(bb, new)

    ob_container = None
    ib_container = None
    for i in range(blanket_final_comp.occurrences.count):
        occ = blanket_final_comp.occurrences.item(i)
        nm = _safe_name(occ).strip().lower()
        cn = ""
        try:
            cn = _safe_name(occ.component).strip().lower()
        except Exception:
            pass
        if nm == "ob" or cn == "ob":
            ob_container = occ
        if nm == "ib" or cn == "ib":
            ib_container = occ

    def breeder_occurrences(container_occ, tag: str):
        out = []
        if not container_occ:
            return out
        for j in range(container_occ.childOccurrences.count):
            o = container_occ.childOccurrences.item(j)
            blob = (_safe_name(o) + " " + _safe_name(o.component)).lower()
            m = re.search(rf"{tag.lower()}\s*([0-9]+)", blob)
            if m:
                out.append((int(m.group(1)), o))
        out.sort(key=lambda t: t[0])
        return out

    for idx, occ in breeder_occurrences(ob_container, "OB"):
        renamed_any = False
        for k in range(occ.childOccurrences.count):
            o2 = occ.childOccurrences.item(k)
            if o2.component.bRepBodies.count > 0:
                rename_bodies_in_comp(o2.component, f"OB{idx}")
                renamed_any = True
        if not renamed_any:
            rename_bodies_in_comp(occ.component, f"OB{idx}")

    for idx, occ in breeder_occurrences(ib_container, "IB"):
        renamed_any = False
        for k in range(occ.childOccurrences.count):
            o2 = occ.childOccurrences.item(k)
            if o2.component.bRepBodies.count > 0:
                rename_bodies_in_comp(o2.component, f"IB{idx}")
                renamed_any = True
        if not renamed_any:
            rename_bodies_in_comp(occ.component, f"IB{idx}")

def rename_bodies_inside_occurrences(root_comp: adsk.fusion.Component, rules):
    all_occs = _collect_occ_tree_from_root(root_comp)
    all_occs.sort(key=lambda o: (o.fullPathName or "").lower())

    def _norm(s: str) -> str:
        s = (s or "").strip().lower()
        s = re.sub(r"\s*\(\d+\)\s*$", "", s)
        return s

    for (hint, prefix) in rules:
        hint_l = (hint or "").strip().lower()
        if not hint_l:
            continue

        matched = [o for o in all_occs if hint_l in _occ_blob(o)]
        if not matched:
            continue

        for idx, occ in enumerate(matched, start=1):
            try:
                comp = occ.component
            except Exception:
                continue

            solids = []
            for i in range(comp.bRepBodies.count):
                b = comp.bRepBodies.item(i)
                try:
                    if b.isSolid:
                        solids.append(b)
                except Exception:
                    pass

            if len(solids) == 1:
                _try_set_name(solids[0], f"{prefix}_{idx}")
            elif len(solids) > 1:
                for bi, b in enumerate(solids, start=1):
                    _try_set_name(b, f"{prefix}_{idx}_{bi}")

def build_visible_solid_inventory(root_comp: adsk.fusion.Component):
    def _count_visible_solids_subtree(occ: adsk.fusion.Occurrence) -> int:
        count = 0
        stack = [occ]
        while stack:
            o = stack.pop()
            try:
                for b in list(o.bRepBodies):
                    try:
                        if b.isSolid and b.isVisible:
                            count += 1
                    except Exception:
                        pass
            except Exception:
                pass
            try:
                for j in range(o.childOccurrences.count):
                    stack.append(o.childOccurrences.item(j))
            except Exception:
                pass
        return count

    inv = {}
    try:
        top = [root_comp.occurrences.item(i) for i in range(root_comp.occurrences.count)]
    except Exception:
        return inv

    top.sort(key=lambda o: (o.fullPathName or "").lower())
    for occ in top:
        try:
            comp_name = occ.component.name
        except Exception:
            comp_name = occ.name if hasattr(occ, "name") else "(unnamed)"
        inv[comp_name or "(unnamed)"] = int(_count_visible_solids_subtree(occ))
    return inv

def write_user_inputs_json(out_path_json: str,
                          root_comp: adsk.fusion.Component,
                          step_export_path: str):
    OB_N_LAYERS = len([x for x in OB_OFFSETS_CM if x < 0]) + 1
    IB_N_LAYERS = len([x for x in IB_OFFSETS_CM if x < 0]) + 1
    N_BREEDER   = len(FRACTIONS) + 1

    inventory = build_visible_solid_inventory(root_comp)

    payload = {
        "meta": {
            "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "units": "cm"
        },
        "geometry": {
            "step_file": step_export_path,
            "armor_cm": float(ARMOR_CM),
            "fw_cm": float(FW_CM),
            "ob_offsets_cm": [float(x) for x in OB_OFFSETS_CM],
            "ib_offsets_cm": [float(x) for x in IB_OFFSETS_CM],
            "fractions": [float(x) for x in FRACTIONS],
            "ob_n_layers": int(OB_N_LAYERS),
            "ib_n_layers": int(IB_N_LAYERS),
            "n_breeder": int(N_BREEDER),
            "ob_regions": 3,
            "ib_regions": 2,
        },
        "inventory": inventory
    }

    with open(out_path_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

# ======================================================
# OB/IB path/plane helpers
# ======================================================
def _edge_param_extents(edge: adsk.fusion.BRepEdge):
    ev = edge.evaluator
    ok, pmin, pmax = ev.getParameterExtents()
    if not ok:
        return 0.0, 1.0
    return pmin, pmax

def _edge_point(edge: adsk.fusion.BRepEdge, t: float) -> adsk.core.Point3D:
    ev = edge.evaluator
    ok, p = ev.getPointAtParameter(t)
    if not ok:
        raise RuntimeError("getPointAtParameter failed")
    return p

def _edge_midpoint(edge: adsk.fusion.BRepEdge) -> adsk.core.Point3D:
    pmin, pmax = _edge_param_extents(edge)
    return _edge_point(edge, 0.5 * (pmin + pmax))

def _edge_length(edge: adsk.fusion.BRepEdge) -> float:
    try:
        return float(edge.length)
    except Exception:
        return 0.0

def _edge_z_span(edge: adsk.fusion.BRepEdge) -> float:
    pmin, pmax = _edge_param_extents(edge)
    zs = []
    for tfrac in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        t = pmin + (pmax - pmin) * tfrac
        try:
            zs.append(_edge_point(edge, t).z)
        except Exception:
            pass
    if len(zs) < 2:
        return 0.0
    return float(max(zs) - min(zs))

def _split_edges_by_radius(edges):
    rs = []
    for e in edges:
        try:
            mp = _edge_midpoint(e)
            rs.append(math.hypot(mp.x, mp.y))
        except Exception:
            rs.append(0.0)
    if not rs:
        return [], [], 0.0
    med = sorted(rs)[len(rs)//2]
    outer, inner = [], []
    for e, r in zip(edges, rs):
        (outer if r >= med else inner).append((r, e))
    return outer, inner, med

def _sample_edge_points(edge: adsk.fusion.BRepEdge, npts: int):
    pmin, pmax = _edge_param_extents(edge)
    if npts < 2:
        npts = 2
    pts = []
    for i in range(npts):
        t = pmin + (pmax - pmin) * (i / (npts - 1))
        pts.append(_edge_point(edge, t))
    if len(pts) < 2:
        raise RuntimeError("Edge sampling produced <2 points.")
    return pts

def _make_root_path_from_3d_spline(root: adsk.fusion.Component, sample_pts_list):
    pts = adsk.core.ObjectCollection.create()
    for p in sample_pts_list:
        pts.add(p)
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.is3D = True
    sp = sk.sketchCurves.sketchFittedSplines.add(pts)
    oc = adsk.core.ObjectCollection.create()
    oc.add(sp)
    path = adsk.fusion.Path.create(oc, adsk.fusion.ChainedCurveOptions.connectedChainedCurves)
    try:
        sk.isVisible = False
    except Exception:
        pass
    return path

def _add_root_plane_along_path(root: adsk.fusion.Component, path: adsk.fusion.Path, frac: float, name: str):
    planes = root.constructionPlanes
    pin = planes.createInput()
    pin.setByDistanceOnPath(path, adsk.core.ValueInput.createByReal(frac))
    pl = planes.add(pin)
    _try_set_name(pl, name)
    return pl

def _delete_existing_planes(root: adsk.fusion.Component, prefix_rx: str):
    planes = root.constructionPlanes
    rx = re.compile(prefix_rx, re.IGNORECASE)
    to_del = []
    for i in range(planes.count):
        p = planes.item(i)
        if rx.match(_safe_name(p)):
            to_del.append(p)
    for p in to_del:
        try:
            p.deleteMe()
        except Exception:
            pass

def _dist2(p: adsk.core.Point3D, q: adsk.core.Point3D) -> float:
    dx = p.x - q.x
    dy = p.y - q.y
    dz = p.z - q.z
    return dx*dx + dy*dy + dz*dz

def _build_path_cumulative_lengths(pts):
    cum = [0.0]
    total = 0.0
    for i in range(1, len(pts)):
        a = pts[i-1]
        b = pts[i]
        seg = math.sqrt((b.x-a.x)**2 + (b.y-a.y)**2 + (b.z-a.z)**2)
        total += seg
        cum.append(total)
    return cum, total

def _centroid_to_path_fraction(com: adsk.core.Point3D, pts, cum, total) -> float:
    best_i = 0
    best_d = 1e99
    for i, p in enumerate(pts):
        d = _dist2(com, p)
        if d < best_d:
            best_d = d
            best_i = i
    if total <= 0.0:
        return 0.0
    return max(0.0, min(1.0, cum[best_i] / total))

def _frac_to_k(frac: float) -> int:
    if frac <= FRACTIONS[0]:
        return 1
    for i in range(1, len(FRACTIONS)):
        if frac <= FRACTIONS[i]:
            return i + 1
    return len(FRACTIONS) + 1

# ======================================================
# IB cuts helpers
# ======================================================
def _collect_solid_bodies_with_occ(occ0: adsk.fusion.Occurrence):
    out = []
    stack = [occ0]
    while stack:
        occ = stack.pop()
        try:
            for b in list(occ.bRepBodies):
                try:
                    if b.isSolid:
                        out.append((b, occ))
                except Exception:
                    pass
        except Exception:
            pass
        try:
            for j in range(occ.childOccurrences.count):
                stack.append(occ.childOccurrences.item(j))
        except Exception:
            pass
    return out

def _brep_edge_length(e: adsk.fusion.BRepEdge) -> float:
    try:
        return float(e.length)
    except Exception:
        return 0.0

def _edges_share_vertex(e1: adsk.fusion.BRepEdge, e2: adsk.fusion.BRepEdge) -> bool:
    try:
        v1a, v1b = e1.startVertex, e1.endVertex
        v2a, v2b = e2.startVertex, e2.endVertex
        return (v1a == v2a) or (v1a == v2b) or (v1b == v2a) or (v1b == v2b)
    except Exception:
        return False

def _order_connected_edges(e1, e2):
    try:
        if e1.endVertex == e2.startVertex:
            return e1, e2
        if e2.endVertex == e1.startVertex:
            return e2, e1
    except Exception:
        pass
    return e1, e2

def _largest_edge_on_body(body: adsk.fusion.BRepBody):
    best = None
    bestL = -1.0
    try:
        for i in range(body.edges.count):
            e = body.edges.item(i)
            L = _brep_edge_length(e)
            if L > bestL:
                best, bestL = e, L
    except Exception:
        pass
    return best, bestL

def _largest_connected_edge(body: adsk.fusion.BRepBody, base_edge: adsk.fusion.BRepEdge):
    best = None
    bestL = -1.0
    try:
        for i in range(body.edges.count):
            e = body.edges.item(i)
            if e == base_edge:
                continue
            if not _edges_share_vertex(base_edge, e):
                continue
            L = _brep_edge_length(e)
            if L > bestL:
                best, bestL = e, L
    except Exception:
        pass
    return best, bestL

def _find_max_layer_index_for_prefix(prefix: str, bodies_with_occ):
    nmax = -1
    rx = re.compile(rf"^{re.escape(prefix)}_Layer_(\d+)$", re.IGNORECASE)
    rx2 = re.compile(rf"^{re.escape(prefix.split('1')[0])}_Layer_(\d+)$", re.IGNORECASE)
    for (b, _) in bodies_with_occ:
        nm = _safe_name(b).strip()
        m = rx.match(nm) or rx2.match(nm)
        if m:
            try:
                nmax = max(nmax, int(m.group(1)))
            except Exception:
                pass
    return nmax

def _find_layer_body_with_occ(prefix: str, n: int, bodies_with_occ):
    rx = re.compile(rf"^{re.escape(prefix)}_Layer_{n}$", re.IGNORECASE)
    rx2 = re.compile(rf"^{re.escape(prefix.split('1')[0])}_Layer_{n}$", re.IGNORECASE)
    for (b, occ) in bodies_with_occ:
        nm = _safe_name(b).strip()
        if rx.match(nm) or rx2.match(nm):
            return b, occ
    return None, None

def _sample_two_connected_edges_points(e1, e2, npts_total: int):
    eA, eB = _order_connected_edges(e1, e2)
    if npts_total < 10:
        npts_total = 10
    nA = max(5, npts_total // 2)
    nB = max(5, npts_total - nA)
    ptsA = _sample_edge_points(eA, nA)
    ptsB = _sample_edge_points(eB, nB)
    if ptsA and ptsB:
        try:
            if _dist2(ptsA[-1], ptsB[0]) < 1e-12:
                ptsB = ptsB[1:]
        except Exception:
            pass
    pts = ptsA + ptsB
    if len(pts) < 2:
        raise RuntimeError("Two-edge sampling produced <2 points.")
    return pts

# ======================================================
# Copy VV pieces into OB_SORTED/IB_SORTED buckets by k
# ======================================================
def copy_vv_into_sorted(vv_occ: adsk.fusion.Occurrence,
                        ob_sorted_comp: adsk.fusion.Component,
                        ib_sorted_comp: adsk.fusion.Component,
                        path_data_by_j,
                        path_data_ib_by_j):
    def _normalize_fusion_name(s: str) -> str:
        s = (s or "").strip().lower()
        s = re.sub(r"\s*\(\d+\)\s*$", "", s)
        return s

    def _find_child_comp_by_name(parent_comp, want_name: str):
        wn = _normalize_fusion_name(want_name)
        best = None
        best_is_exact = False
        for i in range(parent_comp.occurrences.count):
            o = parent_comp.occurrences.item(i)
            try:
                cn = _normalize_fusion_name(_safe_name(o.component))
            except Exception:
                cn = ""
            on = _normalize_fusion_name(_safe_name(o))
            if cn == wn or on == wn:
                raw_on = (_safe_name(o) or "").strip().lower()
                raw_cn = (_safe_name(o.component) or "").strip().lower()
                is_exact = (raw_on == wn) or (raw_cn == wn)
                if best is None or (is_exact and not best_is_exact):
                    best = o.component
                    best_is_exact = is_exact
        return best

    def _find_or_make_k_comp(j_comp: adsk.fusion.Component, k: int):
        kname = f"k{str(k).zfill(2)}"
        kc = _find_child_comp_by_name(j_comp, kname)
        if kc:
            return kc
        _, kc = new_component(j_comp, kname)
        return kc

    def _find_j_comp(sorted_comp: adsk.fusion.Component, j: int, prefix: str):
        want = f"{prefix}{j}".lower()
        jc = _find_child_comp_by_name(sorted_comp, want)
        if jc:
            return jc
        _, jc = new_component(sorted_comp, f"{prefix}{j}")
        return jc

    def _parse_vv_name(name_lower: str):
        m = re.match(r"^vv_(ob|ib)(\d+)(?:_k(\d+))?.*$", name_lower)
        if not m:
            return None, None, None
        side = m.group(1).upper()
        j = int(m.group(2))
        k = int(m.group(3)) if m.group(3) else None
        return side, j, k

    def _k_from_projection(side: str, j: int, body_proxy):
        try:
            com = body_proxy.physicalProperties.centerOfMass
        except Exception:
            return 1

        if side == "OB" and j in path_data_by_j:
            pts, cum, total, _ = path_data_by_j[j]
            try:
                return _frac_to_k(_centroid_to_path_fraction(com, pts, cum, total))
            except Exception:
                return 1

        if side == "IB" and j in path_data_ib_by_j:
            pts, cum, total, _ = path_data_ib_by_j[j]
            try:
                return _frac_to_k(_centroid_to_path_fraction(com, pts, cum, total))
            except Exception:
                return 1

        return 1

    vv_bodies = list(_solid_proxy_bodies(vv_occ))
    if not vv_bodies:
        return 0, 0, 0

    copied = 0
    hidden = 0
    skipped = 0

    for b in vv_bodies:
        nm = _safe_name(b).strip()
        nml = nm.lower()

        if nml.startswith("vv_port_") or nml == "vv_divertor":
            continue

        side, j, k = _parse_vv_name(nml)
        if side is None:
            continue
        if side == "OB" and j not in (1, 2, 3):
            continue
        if side == "IB" and j not in (1, 2):
            continue

        if k is None:
            k = _k_from_projection(side, j, b)

        if side == "OB":
            j_comp = _find_j_comp(ob_sorted_comp, j, "OB")
            k_comp = _find_or_make_k_comp(j_comp, k)
        else:
            j_comp = _find_j_comp(ib_sorted_comp, j, "IB")
            k_comp = _find_or_make_k_comp(j_comp, k)

        new_name = f"VV_{side}{j}_k{str(k).zfill(2)}"
        try:
            copy_body_into_def(k_comp, b, new_name)
            copied += 1
            try:
                b.isVisible = False
                hidden += 1
            except Exception:
                pass
        except Exception:
            skipped += 1

    return copied, hidden, skipped

# ======================================================
# Port filling (FIXED width/height swap)
# ======================================================
def create_port_filling_block(root_comp: adsk.fusion.Component,
                              vv_occ: adsk.fusion.Occurrence,
                              port_body_name: str = "VV_port_2",
                              theta_deg: float = 11.25,
                              width_cm: float = 180.0,
                              height_cm: float = 250.0,
                              swap_wh: bool = False,
                              use_square: bool = False,
                              depth_scale: float = 0.60):
    try:
        port_body = _find_body_by_exact_name_in_occ(vv_occ, port_body_name)
        if port_body is None:
            return False, f"port_filling: FAIL - '{port_body_name}' not found.", None

        com = port_body.physicalProperties.centerOfMass

        th = math.radians(theta_deg)
        ax = math.cos(th)
        ay = math.sin(th)
        az = 0.0

        pf_occ = root_comp.occurrences.addNewComponent(adsk.core.Matrix3D.create())
        pf_comp = pf_occ.component
        _try_set_name(pf_comp, "PORT_FILLING")

        planes_root = root_comp.constructionPlanes

        pin_rot = planes_root.createInput()
        pin_rot.setByAngle(root_comp.zConstructionAxis,
                           adsk.core.ValueInput.createByReal(th),
                           root_comp.yZConstructionPlane)
        pl_rot = planes_root.add(pin_rot)
        _try_set_name(pl_rot, "PORT_FILLING_PLANE_ROT")

        dist = com.x * ax + com.y * ay + com.z * az
        pin_off = planes_root.createInput()
        pin_off.setByOffset(pl_rot, adsk.core.ValueInput.createByReal(dist))
        pl_at_com = planes_root.add(pin_off)
        _try_set_name(pl_at_com, "PORT_FILLING_PLANE_AT_COM")

        try:
            pl_pf = pl_at_com.createForAssemblyContext(pf_occ)
        except Exception:
            pl_pf = pl_at_com

        sk = pf_comp.sketches.add(pl_pf)
        c2d = sk.modelToSketchSpace(com)

        w = float(width_cm)
        h = float(height_cm if not use_square else width_cm)

        # --- FIX: if swapped in your Fusion view, swap drawing axes ---
        if swap_wh:
            w, h = h, w

        if hasattr(sk.sketchCurves, "sketchLines") and hasattr(sk.sketchCurves.sketchLines, "addCenterPointRectangle"):
            sk.sketchCurves.sketchLines.addCenterPointRectangle(
                adsk.core.Point3D.create(c2d.x, c2d.y, 0),
                adsk.core.Point3D.create(c2d.x + 0.5*w, c2d.y + 0.5*h, 0)
            )
        elif hasattr(sk, "sketchLines") and hasattr(sk.sketchLines, "addCenterPointRectangle"):
            sk.sketchLines.addCenterPointRectangle(
                adsk.core.Point3D.create(c2d.x, c2d.y, 0),
                adsk.core.Point3D.create(c2d.x + 0.5*w, c2d.y + 0.5*h, 0)
            )
        else:
            return False, "port_filling: FAIL - no rectangle API available in this Fusion build.", None

        if sk.profiles.count < 1:
            return False, "port_filling: FAIL - sketch produced no profiles.", None
        prof = sk.profiles.item(0)

        # depth estimate by bbox projection on axis
        bb = port_body.boundingBox
        corners = [
            adsk.core.Point3D.create(bb.minPoint.x, bb.minPoint.y, bb.minPoint.z),
            adsk.core.Point3D.create(bb.minPoint.x, bb.minPoint.y, bb.maxPoint.z),
            adsk.core.Point3D.create(bb.minPoint.x, bb.maxPoint.y, bb.minPoint.z),
            adsk.core.Point3D.create(bb.minPoint.x, bb.maxPoint.y, bb.maxPoint.z),
            adsk.core.Point3D.create(bb.maxPoint.x, bb.minPoint.y, bb.minPoint.z),
            adsk.core.Point3D.create(bb.maxPoint.x, bb.minPoint.y, bb.maxPoint.z),
            adsk.core.Point3D.create(bb.maxPoint.x, bb.maxPoint.y, bb.minPoint.z),
            adsk.core.Point3D.create(bb.maxPoint.x, bb.maxPoint.y, bb.maxPoint.z),
        ]
        def dotp(p):
            return p.x*ax + p.y*ay + p.z*az
        projs = [dotp(p) for p in corners]
        Laxis = max(projs) - min(projs)
        if Laxis <= 1e-6:
            Laxis = 200.0
        depth = max(20.0, float(depth_scale) * Laxis)

        ext_feats = pf_comp.features.extrudeFeatures
        ext_in = ext_feats.createInput(prof, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)

        if hasattr(ext_in, "setSymmetricExtent"):
            ext_in.setSymmetricExtent(adsk.core.ValueInput.createByReal(0.5*depth), True)
        else:
            ext_in.setDistanceExtent(False, adsk.core.ValueInput.createByReal(depth))

        ext = ext_feats.add(ext_in)
        if not ext or ext.bodies.count < 1:
            return False, "port_filling: FAIL - extrude produced no bodies.", None

        block = ext.bodies.item(0)
        _try_set_name(block, f"{port_body_name}_FILL")

        return True, (
            f"port_filling: OK - created block, requested={width_cm:.0f}x{height_cm:.0f} cm "
            f"(drawn={w:.0f}x{h:.0f}), depth={depth:.1f} cm, theta={theta_deg}°, swap_wh={swap_wh}"
        ), block

    except Exception as e:
        return False, f"port_filling: FAIL - exception: {e}", None

# ======================================================
# PFC merge triplets + TF join (your earlier requests)
# ======================================================
def merge_pfc_in_place_triplets(root_comp: adsk.fusion.Component,
                               start_idx: int = 19,
                               end_idx: int = 33,
                               group_size: int = 3,
                               keep_tools: bool = False):
    def _norm(s: str) -> str:
        s = (s or "").strip().lower()
        s = re.sub(r"\s*\(\d+\)\s*$", "", s)
        return s

    all_occs = _collect_occ_tree_from_root(root_comp)
    name_to_bodies = {}
    for occ in all_occs:
        for b in list(_solid_proxy_bodies(occ)):
            nm = _norm(_safe_name(b))
            if nm.startswith("pfc_1_"):
                name_to_bodies.setdefault(nm, []).append(b)

    cmb = root_comp.features.combineFeatures

    merged_groups = 0
    join_ops = 0
    missing = []
    fails = 0
    first_err = None

    idx = start_idx
    while idx <= end_idx:
        grp = list(range(idx, min(idx + group_size, end_idx + 1)))
        if len(grp) < group_size:
            break

        want = [f"pfc_1_{k}" for k in grp]
        bodies = []
        for wn in want:
            cands = name_to_bodies.get(wn, [])
            if not cands:
                missing.append(wn)
                continue
            bodies.append(sorted(cands, key=_proxy_volume, reverse=True)[0])

        if len(bodies) != group_size:
            idx += group_size
            continue

        target = sorted(bodies, key=_proxy_volume, reverse=True)[0]
        tools_list = [b for b in bodies if b != target]

        tools = adsk.core.ObjectCollection.create()
        for t in tools_list:
            tools.add(t)

        group_name = f"PFC_1_{grp[0]}_{grp[-1]}"
        try:
            ci = cmb.createInput(target, tools)
            ci.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
            ci.isKeepToolBodies = bool(keep_tools)
            cmb.add(ci)
            _rename_body_native_if_possible(target, group_name)
            merged_groups += 1
            join_ops += 1
        except Exception as e:
            fails += 1
            if first_err is None:
                first_err = str(e)

        idx += group_size

    msg = f"PFC triplet JOIN: merged_groups={merged_groups}, fails={fails}, missing={len(missing)}"
    if first_err:
        msg += f". First error: {first_err}"
    return (fails == 0), msg

def join_all_tfcoils(root_comp: adsk.fusion.Component,
                     tf_hint: str = "tfcoil_mat_toroidal_field_coil_2015",
                     keep_tools: bool = False):
    occs = _find_all_occurrences_by_hint(root_comp, tf_hint)
    if not occs:
        return False, f"TFcoil join: no occurrences found for hint '{tf_hint}'."

    solids = []
    for occ in occs:
        solids.extend(list(_solid_proxy_bodies(occ)))

    if len(solids) < 2:
        return True, f"TFcoil join: found {len(solids)} solid(s); nothing to join."

    solids_sorted = sorted(solids, key=_proxy_volume, reverse=True)
    target = solids_sorted[0]
    tools_list = solids_sorted[1:]

    tools = adsk.core.ObjectCollection.create()
    for b in tools_list:
        tools.add(b)

    try:
        cmb = root_comp.features.combineFeatures
        ci = cmb.createInput(target, tools)
        ci.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
        ci.isKeepToolBodies = bool(keep_tools)
        cmb.add(ci)
        _rename_body_native_if_possible(target, "TFcoil_merged")
        return True, f"TFcoil join: joined {len(tools_list)} body(ies) into TFcoil_merged. keep_tools={keep_tools}"
    except Exception as e:
        return False, f"TFcoil join: JOIN failed: {e}"

# ======================================================
# Move VV_ports
# ======================================================
def move_vv_ports(root_comp: adsk.fusion.Component,
                  vv_occ: adsk.fusion.Occurrence,
                  dx_cm: float = 2.0,
                  dz_cm: float = 5.0,
                  port_prefix: str = "VV_port_"):
    try:
        if vv_occ is None:
            return (False, "move_vv_ports: vv_occ is None")

        port_bodies = []
        for b in list(_solid_proxy_bodies(vv_occ)):
            nm = _safe_name(b).strip()
            if nm.lower().startswith(port_prefix.lower()):
                port_bodies.append(b)

        if not port_bodies:
            return (True, "move_vv_ports: no VV_port_* bodies found; nothing to move.")

        def _desired_delta_for_name(nm: str):
            lo = nm.strip().lower()
            if lo.startswith("vv_port_3"):
                return (0.0, 0.0, float(dz_cm))
            return (float(dx_cm), 0.0, 0.0)

        def _to_native_body(proxy_body: adsk.fusion.BRepBody):
            try:
                return proxy_body.nativeObject if proxy_body and proxy_body.isValid else None
            except Exception:
                return None

        moves_by_owner = {}

        for pb in port_bodies:
            nm = _safe_name(pb)
            nb = _to_native_body(pb) or pb
            try:
                owner_comp = nb.parentComponent
            except Exception:
                owner_comp = None

            if owner_comp is not None:
                try:
                    key = owner_comp.entityToken
                except Exception:
                    key = f"{_safe_name(owner_comp)}::{id(owner_comp)}"
            else:
                key = f"UNKNOWN::{id(nb)}"

            dxyz = _desired_delta_for_name(nm)
            moves_by_owner.setdefault(key, {"comp": owner_comp, "items": []})["items"].append((nb, dxyz, nm))

        ok_count = 0
        fail_count = 0
        lines = []

        for _, entry in moves_by_owner.items():
            owner_comp = entry["comp"]
            items = entry["items"]

            if owner_comp is None or (hasattr(owner_comp, "isValid") and not owner_comp.isValid):
                for (_, _, nm) in items:
                    fail_count += 1
                    lines.append(f"  '{nm}': FAIL (no valid owner component)")
                continue

            mf = owner_comp.features.moveFeatures

            for (nb, dxyz, nm) in items:
                try:
                    dx, dy, dz = dxyz
                    vec = adsk.core.Vector3D.create(dx, dy, dz)
                    tr = adsk.core.Matrix3D.create()
                    tr.translation = vec

                    objs = adsk.core.ObjectCollection.create()
                    objs.add(nb)

                    inp = mf.createInput(objs, tr)
                    mf.add(inp)

                    ok_count += 1
                    lines.append(f"  '{nm}': moved approx ({dx:+.1f},{dy:+.1f},{dz:+.1f}) cm")
                except Exception as e_one:
                    fail_count += 1
                    lines.append(f"  '{nm}': FAIL ({e_one})")

        ok_all = (fail_count == 0)
        msg = "move_vv_ports: " + ("OK" if ok_all else "FAIL") + f" ok={ok_count}, fail={fail_count}\n" + "\n".join(lines)
        return (ok_all, msg)

    except Exception as e:
        return (False, f"move_vv_ports: exception: {e}")

# ======================================================
# 15 de-overlap, BF x VV CUT, and overlap verify
# ======================================================
def safe_body_volume_cm3(body):
    if body is None:
        return (False, 0.0, "body is None")
    try:
        v = float(body.physicalProperties.volume)
        return (v > 0.0), v, "ok" if v > 0 else "volume=0"
    except Exception as e:
        return (False, 0.0, f"physicalProperties failed: {e}")

def _safe_is_visible(obj, default=True):
    try:
        return bool(obj.isVisible)
    except Exception:
        return default

def _safe_is_solid(body):
    try:
        return bool(body.isSolid)
    except Exception:
        return True

def _body_volume(body) -> float:
    try:
        return float(body.physicalProperties.volume)
    except Exception:
        return 0.0

def _collect_visible_solids_in_comp_def(comp: adsk.fusion.Component):
    out = []
    for i in range(comp.bRepBodies.count):
        b = comp.bRepBodies.item(i)
        if not _safe_is_solid(b):
            continue
        if not _safe_is_visible(b, default=True):
            continue
        out.append(b)
    return out

def _combine_cut_keep_tools(comp: adsk.fusion.Component,
                           target: adsk.fusion.BRepBody,
                           tools_list,
                           tools_chunk: int = 40):
    cmb = comp.features.combineFeatures
    ok = 0
    fail = 0

    filtered = []
    for t in tools_list:
        if t is None or t == target:
            continue
        if not _safe_is_solid(t):
            continue
        try:
            if hasattr(t, "isVisible") and not t.isVisible:
                continue
        except Exception:
            continue
        filtered.append(t)

    if not filtered:
        return 0, 0

    for i0 in range(0, len(filtered), max(1, int(tools_chunk))):
        chunk = filtered[i0:i0 + tools_chunk]
        oc = adsk.core.ObjectCollection.create()
        for t in chunk:
            oc.add(t)
        try:
            ci = cmb.createInput(target, oc)
            ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
            ci.isKeepToolBodies = True
            cmb.add(ci)
            ok += 1
        except Exception:
            fail += 1
    return ok, fail

def boolean_deoverlap_blanket_final(bf_comp: adsk.fusion.Component,
                                    min_keep_vol_cm3: float = 1e-6,
                                    tools_chunk: int = 40):
    bodies = _collect_visible_solids_in_comp_def(bf_comp)
    if len(bodies) < 2:
        return True, f"13.b deoverlap: <2 visible solids in '{_safe_name(bf_comp)}'.", {}

    bodies_sorted = sorted(bodies, key=_body_volume, reverse=True)
    accepted = []
    cut_ok = 0
    cut_fail = 0
    deleted_zero = 0

    for b in bodies_sorted:
        if not _safe_is_solid(b) or not _safe_is_visible(b, default=True):
            continue

        if accepted:
            ok, fail = _combine_cut_keep_tools(bf_comp, b, accepted, tools_chunk=tools_chunk)
            cut_ok += ok
            cut_fail += fail

        try:
            v = _body_volume(b)
            if v < float(min_keep_vol_cm3):
                try:
                    b.deleteMe()
                except Exception:
                    try:
                        b.isVisible = False
                    except Exception:
                        pass
                deleted_zero += 1
                continue
        except Exception:
            pass

        accepted.append(b)

    msg = (f"13.b deoverlap: processed {len(bodies)} visible solids; "
           f"accepted={len(accepted)}, cut_ops_ok={cut_ok}, cut_ops_fail={cut_fail}, deleted_zero={deleted_zero}.")
    stats = {"visible_solids_start": len(bodies), "accepted_end": len(accepted),
             "cut_ops_ok": cut_ok, "cut_ops_fail": cut_fail, "deleted_zero": deleted_zero}
    return True, msg, stats

def _collect_visible_proxy_solids_under_occ(occ0: adsk.fusion.Occurrence):
    out = []
    stack = [occ0]
    while stack:
        occ = stack.pop()
        try:
            for b in list(occ.bRepBodies):
                try:
                    if not b.isSolid:
                        continue
                    if hasattr(b, "isVisible") and not b.isVisible:
                        continue
                    out.append(b)
                except Exception:
                    pass
        except Exception:
            pass
        try:
            for j in range(occ.childOccurrences.count):
                stack.append(occ.childOccurrences.item(j))
        except Exception:
            pass
    return out

def _combine_cut_keep_tools_root(root_comp: adsk.fusion.Component,
                                target_proxy: adsk.fusion.BRepBody,
                                tools_proxies,
                                tools_chunk: int = 40):
    cmb = root_comp.features.combineFeatures
    ok = 0
    fail = 0

    filtered = []
    for t in tools_proxies:
        if t is None or t == target_proxy:
            continue
        try:
            if not t.isSolid:
                continue
        except Exception:
            pass
        try:
            if hasattr(t, "isVisible") and not t.isVisible:
                continue
        except Exception:
            continue
        filtered.append(t)

    if not filtered:
        return 0, 0

    for i0 in range(0, len(filtered), max(1, int(tools_chunk))):
        chunk = filtered[i0:i0 + tools_chunk]
        oc = adsk.core.ObjectCollection.create()
        for t in chunk:
            oc.add(t)
        try:
            ci = cmb.createInput(target_proxy, oc)
            ci.operation = adsk.fusion.FeatureOperations.CutFeatureOperation
            ci.isKeepToolBodies = True
            cmb.add(ci)
            ok += 1
        except Exception:
            fail += 1
    return ok, fail

def boolean_cut_blanket_final_by_vacuum_vessel_visible(root_comp: adsk.fusion.Component,
                                                       bf_occ: adsk.fusion.Occurrence,
                                                       vv_occ: adsk.fusion.Occurrence,
                                                       tools_chunk: int = 40,
                                                       min_keep_vol_cm3: float = 1e-6,
                                                       vv_name_exclude_prefixes=("vv_port_",),
                                                       vv_name_exclude_exact=("vv_divertor",)):
    if bf_occ is None or vv_occ is None:
        return False, "BF x VV CUT: missing occurrence(s).", {}

    bf_targets = _collect_visible_proxy_solids_under_occ(bf_occ)
    if not bf_targets:
        return True, "BF x VV CUT: no visible BLANKET_FINAL solids.", {"bf_targets": 0, "vv_tools": 0}

    vv_tools_all = _collect_visible_proxy_solids_under_occ(vv_occ)
    vv_tools = []
    for b in vv_tools_all:
        nm = (_safe_name(b) or "").strip().lower()
        if any(nm.startswith(p) for p in (vv_name_exclude_prefixes or ())):
            continue
        if nm in set((vv_name_exclude_exact or ())):
            continue
        vv_tools.append(b)

    if not vv_tools:
        return True, "BF x VV CUT: no visible VV tool solids (after excludes).", {"bf_targets": len(bf_targets), "vv_tools": 0}

    cut_ok = 0
    cut_fail = 0
    deleted_zero = 0

    bf_targets_sorted = sorted(bf_targets, key=_proxy_volume, reverse=True)
    for tgt in bf_targets_sorted:
        try:
            if hasattr(tgt, "isVisible") and not tgt.isVisible:
                continue
        except Exception:
            continue

        ok, fail = _combine_cut_keep_tools_root(root_comp, tgt, vv_tools, tools_chunk=tools_chunk)
        cut_ok += ok
        cut_fail += fail

        try:
            okv, v, _ = safe_body_volume_cm3(tgt)
            if okv and v < float(min_keep_vol_cm3):
                try:
                    tgt.deleteMe()
                except Exception:
                    try:
                        tgt.isVisible = False
                    except Exception:
                        pass
                deleted_zero += 1
        except Exception:
            pass

    stats = {"bf_targets": len(bf_targets), "vv_tools": len(vv_tools),
             "cut_ops_ok": cut_ok, "cut_ops_fail": cut_fail, "deleted_zero": deleted_zero}
    msg = (f"BF x VV CUT: targets={len(bf_targets)} tools={len(vv_tools)} "
           f"cut_ops_ok={cut_ok}, cut_ops_fail={cut_fail}, deleted_zero={deleted_zero}.")
    return True, msg, stats

def detect_overlaps_in_component(comp: adsk.fusion.Component,
                                 min_overlap_cm3: float = 1e-6,
                                 only_visible: bool = True,
                                 max_pairs: int = 20000):
    def _solids(c):
        out = []
        for i in range(c.bRepBodies.count):
            b = c.bRepBodies.item(i)
            try:
                if not b.isSolid:
                    continue
                if only_visible and hasattr(b, "isVisible") and not b.isVisible:
                    continue
                out.append(b)
            except Exception:
                pass
        return out

    bodies = _solids(comp)
    n = len(bodies)
    if n < 2:
        return [], f"overlap_check: '{_safe_name(comp)}' has <2 solids."

    pairs = n * (n - 1) // 2
    if pairs > max_pairs:
        return [], f"overlap_check: too many pairs ({pairs}); set max_pairs or filter."

    cmb = comp.features.combineFeatures
    cps = comp.features.copyPasteBodies

    overlaps = []
    checked = 0
    for i in range(n):
        a = bodies[i]
        for j in range(i + 1, n):
            b = bodies[j]
            checked += 1

            fa = cps.add(a)
            if not fa or fa.bodies.count < 1:
                continue
            a2 = fa.bodies.item(0)

            tools = adsk.core.ObjectCollection.create()
            tools.add(b)
            try:
                ci = cmb.createInput(a2, tools)
                ci.operation = adsk.fusion.FeatureOperations.IntersectFeatureOperation
                ci.isKeepToolBodies = True
                cmb.add(ci)
            except Exception:
                try:
                    a2.deleteMe()
                except Exception:
                    pass
                continue

            okv, v, _ = safe_body_volume_cm3(a2)
            if okv and v >= float(min_overlap_cm3):
                overlaps.append((_safe_name(a), _safe_name(b), float(v)))

            try:
                a2.deleteMe()
            except Exception:
                pass

    overlaps.sort(key=lambda t: t[2], reverse=True)
    if overlaps:
        return overlaps, f"overlap_check: FOUND {len(overlaps)} overlaps (checked {checked} pairs)."
    return [], f"overlap_check: none above {min_overlap_cm3:g} (checked {checked} pairs)."

# ======================================================
# Main
# ======================================================
def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        step_timer = StepTimer()
        _log_to_text_palette(app, "Tokamak_test_2 (clean) started...\n")

        if not os.path.isfile(STEP_PATH):
            ui.messageBox(f"STEP file not found:\n{STEP_PATH}")
            return

        # 1) Import STEP
        import_mgr = app.importManager
        step_opts = import_mgr.createSTEPImportOptions(STEP_PATH)
        step_opts.isViewFit = False

        doc = import_mgr.importToNewDocument(step_opts)
        if not doc:
            ui.messageBox("importToNewDocument returned None (import failed).")
            return

        doc.activate()
        design = adsk.fusion.Design.cast(app.activeProduct)

        global root
        root = design.rootComponent
        step_timer.mark("01 Import STEP")

        # 2) Units
        if SET_DISPLAY_UNITS_TO_CM:
            try:
                design.unitsManager.distanceDisplayUnits = adsk.fusion.DistanceUnits.CentimeterDistanceUnits
            except Exception:
                pass

        # 3) Plasma
        def create_plasma_region_torus_pipe(root_comp, major_radius_cm, tube_radius_cm):
            occ = root_comp.occurrences.addNewComponent(adsk.core.Matrix3D.create())
            comp = occ.component
            _try_set_name(comp, "Plasma_Region")

            sk_path = comp.sketches.add(comp.xYConstructionPlane)
            path_circle = sk_path.sketchCurves.sketchCircles.addByCenterRadius(
                adsk.core.Point3D.create(0, 0, 0), major_radius_cm
            )
            path = comp.features.createPath(path_circle)

            pipe_feats = comp.features.pipeFeatures
            pipe_in = pipe_feats.createInput(path, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
            if hasattr(pipe_in, "isSolid"):
                pipe_in.isSolid = True
            if hasattr(pipe_in, "sectionSize"):
                pipe_in.sectionSize = adsk.core.ValueInput.createByReal(2.0 * tube_radius_cm)
            elif hasattr(pipe_in, "radius"):
                pipe_in.radius = adsk.core.ValueInput.createByReal(tube_radius_cm)
            elif hasattr(pipe_in, "diameter"):
                pipe_in.diameter = adsk.core.ValueInput.createByReal(2.0 * tube_radius_cm)
            else:
                raise RuntimeError("PipeInput missing sectionSize/radius/diameter")

            pipe_feat = pipe_feats.add(pipe_in)
            body = None
            if pipe_feat and pipe_feat.bodies.count > 0:
                body = pipe_feat.bodies.item(0)
                _try_set_name(body, "Plasma_Region")
            return occ, comp, body

        #plasma_occ, plasma_comp, plasma_body = create_plasma_region_torus_pipe(
        #    root, PLASMA_MAJOR_RADIUS_CM, PLASMA_TUBE_RADIUS_CM
        #)

        plasma_occ = plasma_comp = plasma_body = None

        if DO_CREATE_PLASMA_REGION:
            plasma_occ, plasma_comp, plasma_body = create_plasma_region_torus_pipe(
                root, PLASMA_MAJOR_RADIUS_CM, PLASMA_TUBE_RADIUS_CM
            )
            step_timer.mark("03 Create Plasma_Region")
        else:
            _log_to_text_palette(app, "03 Create Plasma_Region: skipped (DO_CREATE_PLASMA_REGION=False)\n")
            step_timer.mark("03 Create Plasma_Region (skipped)")

        step_timer.mark("03 Create Plasma_Region")

        # 4) Wedge planes
        plane0 = root.xZConstructionPlane
        planes = root.constructionPlanes
        p_in = planes.createInput()
        p_in.setByAngle(root.zConstructionAxis,
                        adsk.core.ValueInput.createByReal(math.radians(ANGLE_MAX_DEG)),
                        plane0)
        planeA = planes.add(p_in)
        _try_set_name(planeA, f"WedgePlane_{ANGLE_MAX_DEG}deg_ROOT")
        step_timer.mark("04 Create wedge planes")

        # 5) Split by wedge planes
        split_feats_root = root.features.splitBodyFeatures

        ok0 = sk0 = 0
        occs = _collect_occ_tree_from_root(root)
        for occ in occs:
            for body in list(_solid_proxy_bodies(occ)):
                try:
                    split_feats_root.add(split_feats_root.createInput(body, plane0, True))
                    ok0 += 1
                except Exception:
                    sk0 += 1

        okA = skA = 0
        occs = _collect_occ_tree_from_root(root)  # refresh after plane0 split
        for occ in occs:
            for body in list(_solid_proxy_bodies(occ)):
                try:
                    split_feats_root.add(split_feats_root.createInput(body, planeA, True))
                    okA += 1
                except Exception:
                    skA += 1

        step_timer.mark("05 Split by wedge planes")

        # 6) Remove outside wedge
        remove_feats = root.features.removeFeatures
        removed = removed_low = removed_high = remove_fail = 0

        occs = _collect_occ_tree_from_root(root)
        for occ in occs:
            for b in _solid_proxy_bodies(occ):
                try:
                    com = b.physicalProperties.centerOfMass
                except Exception:
                    continue
                theta = math.degrees(math.atan2(com.y, com.x))
                if theta < (0.0 - ANG_TOL_DEG) or theta > (ANGLE_MAX_DEG + ANG_TOL_DEG):
                    try:
                        remove_feats.add(b)
                        removed += 1
                        if theta < 0:
                            removed_low += 1
                        else:
                            removed_high += 1
                    except Exception:
                        remove_fail += 1

        step_timer.mark("06 Remove outside wedge")

        # 7) Delete named occurrences
        deleted_occ = delete_fail = 0
        occs = _collect_occ_tree_from_root(root)
        for occ in sorted(occs, key=lambda o: o.fullPathName.count("+"), reverse=True):
            if not _name_matches_remove_list(occ):
                continue
            try:
                occ.deleteMe()
                deleted_occ += 1
            except Exception:
                delete_fail += 1

        step_timer.mark("07 Delete named occurrences")

        # 8) Rename non-breeder components bodies
        rename_bodies_inside_occurrences(root, RENAME_RULES)

        # 9a) Merge PFC triplets
        if DO_PFC_TRIPLET_MERGE:
            pfc_ok, pfc_msg = merge_pfc_in_place_triplets(root)
            _log_to_text_palette(app, pfc_msg + "\n")

        # 9b) Merge TFcoils
        if DO_TFCOIL_JOIN:
            tf_ok, tf_msg = join_all_tfcoils(root)
            _log_to_text_palette(app, tf_msg + "\n")

        # 10) VV operations
        vv_report = []
        vv_occ = _find_first_occurrence_by_hint(root, VV_HINT)
        ob_occ_pre = _find_first_occurrence_by_hint(root, OB_HINT)
        ib_occ_pre = _find_first_occurrence_by_hint(root, IB_HINT)

        if DO_VV_OPS_EARLY and vv_occ and ob_occ_pre and ib_occ_pre:
            # 10.1) Split VV into VV_1 + VV_port_n (or just rename if already split)
            _, _, vv_msg = split_and_name_vv_ports(root, VV_HINT)
            vv_report.append(vv_msg)

            # refresh occurrence handle (safe)
            vv_occ = _find_first_occurrence_by_hint(root, VV_HINT)

            # 10.2) Move ports (your newer preference)
            mv_ok, mv_msg = move_vv_ports(root, vv_occ, dx_cm=2.0, dz_cm=5.0)
            vv_report.append(mv_msg)

            # 10.3) Face cleanup on VV_1 (keep largest N)
            ok_clean, msg_clean = vv1_face_cleanup_keep_largest(
                vv_occ, keep_faces=VV_KEEP_FACES, min_face_area=VV_MIN_FACE_AREA, heal=VV_HEAL_ON_DELETE
            )
            vv_report.append(msg_clean)

            # 10.4) CUT ports against VV_1 (ports := ports - VV_1) -> Most likely this could be removed
            ok_cut, msg_cut = vv_ports_cut_against_vv1(root, vv_occ, keep_tool_vv1=True)
            vv_report.append(msg_cut)

            # 10.5) Eun these AFTER VV_1 exists and ports are handled.
            vv_occ = _find_first_occurrence_by_hint(root, VV_HINT)
            if vv_occ:
                # 10.5a) Poloidal split directed => VV_OB, VV_IB, VV_divertor
                vvpol_ok, vvpol_msg = vv_poloidal_split_directed(root, vv_occ, ob_occ_pre, ib_occ_pre)
                vv_report.append(vvpol_msg)

                # 10.5b) Split VV_OB into 3 => VV_OB1,VV_OB2,VV_OB3
                ob_ok, ob_msg = vv_split_ob_into_3(root, vv_occ, ob_occ_pre)
                vv_report.append(ob_msg)

                # 10.5c) Split VV_IB into 2 => VV_IB1,VV_IB2
                ib1_occ = ib_occ_pre
                ib2_occ = ib_occ_pre
                ib_ok, ib_msg = vv_split_ib_into_2_by_facing_face(root, vv_occ, ib1_occ, ib2_occ)
                vv_report.append(ib_msg)

        else:
            vv_report.append("VV ops: skipped (VV/OB/IB not found or toggle off).")

        step_timer.mark("08 Early VV ops (FULL)")

        # 11) Blanket layering
        if ob_occ_pre is None:
            raise RuntimeError(f"Could not find OB occurrence matching hint: '{OB_HINT}'")
        if ib_occ_pre is None:
            raise RuntimeError(f"Could not find IB occurrence matching hint: '{IB_HINT}'")

        ob_bodies = _solid_definition_bodies(ob_occ_pre)
        ib_bodies = _solid_definition_bodies(ib_occ_pre)
        if not ob_bodies:
            raise RuntimeError("No solid OB bodies found.")
        if not ib_bodies:
            raise RuntimeError("No solid IB bodies found.")

        blanket_final_occ, blanket_final_comp = new_component(root, "BLANKET_FINAL")
        ob_summary, _ = build_layers_for_breeder(root, blanket_final_comp, ob_bodies, "OB", OB_NTH_FACE, OB_OFFSETS_CM)
        ib_summary, _ = build_layers_for_breeder(root, blanket_final_comp, ib_bodies, "IB", IB_NTH_FACE, IB_OFFSETS_CM)

        if DELETE_ORIGINAL_BLANKET_OCCURRENCES:
            to_delete = []
            to_delete += _find_all_occurrences_by_hint(root, OB_HINT)
            to_delete += _find_all_occurrences_by_hint(root, IB_HINT)
            for occ in sorted(to_delete, key=lambda o: o.fullPathName.count("+"), reverse=True):
                try:
                    occ.deleteMe()
                except Exception:
                    try:
                        occ.isLightBulbOn = False
                    except Exception:
                        pass

        rename_blanket_final(blanket_final_comp)

        blanket_final_occ2 = _find_first_occurrence_by_hint(root, FINAL_COMP_NAME_HINT)
        if blanket_final_occ2 is None:
            raise RuntimeError(f"Could not find occurrence containing '{FINAL_COMP_NAME_HINT}'")

        step_timer.mark("09 Build OB/IB layers")

        # 12.1) OB PLANES (MASTER=OB2)
        ob_occ_root = None
        for i in range(blanket_final_occ2.childOccurrences.count):
            c = blanket_final_occ2.childOccurrences.item(i)
            if _safe_name(c.component).strip().lower() == "ob" or _safe_name(c).strip().lower() == "ob":
                ob_occ_root = c
                break
        if ob_occ_root is None:
            raise RuntimeError("Could not find 'OB' occurrence under BLANKET_FINAL.")

        obj_occs = []
        for i in range(ob_occ_root.childOccurrences.count):
            o = ob_occ_root.childOccurrences.item(i)
            blob = (_safe_name(o) + " " + _safe_name(o.component)).lower()
            m = re.search(r"ob\s*([0-9]+)", blob)
            if m:
                obj_occs.append((int(m.group(1)), o))
        obj_occs.sort(key=lambda t: t[0])
        if not obj_occs:
            raise RuntimeError("No OB{i} occurrences found under BLANKET_FINAL/OB.")

        def _body_com_radius(body):
            com = body.physicalProperties.centerOfMass
            return math.hypot(com.x, com.y)

        _delete_existing_planes(root, r"^OB\d+_PLANE_\d{4}$")

        MASTER_OB_J = 2
        ob_occ_master = None
        for (j, occj) in obj_occs:
            if j == MASTER_OB_J:
                ob_occ_master = occj
                break
        if ob_occ_master is None:
            ob_occ_master = obj_occs[0][1]
            MASTER_OB_J = obj_occs[0][0]

        master_prefix = f"OB{MASTER_OB_J}"
        nmax = -1
        for bb in list(ob_occ_master.bRepBodies):
            try:
                if not bb.isSolid:
                    continue
            except Exception:
                continue
            nm = _safe_name(bb).strip().lower()
            m1 = re.match(rf"^{master_prefix.lower()}_layer_(\d+)$", nm)
            if m1:
                nmax = max(nmax, int(m1.group(1)))
        if nmax < 0:
            raise RuntimeError(f"OB planes: MASTER {master_prefix} has no Layer_N bodies.")

        candidates = []
        for bb in list(ob_occ_master.bRepBodies):
            try:
                if not bb.isSolid:
                    continue
            except Exception:
                continue
            nm = _safe_name(bb).strip().lower()
            if nm == f"{master_prefix.lower()}_layer_{nmax}":
                candidates.append(bb)

        if not candidates:
            raise RuntimeError(f"OB planes: MASTER {master_prefix} missing Layer_{nmax} body.")

        ref_body = None
        best_r = -1.0
        for b in candidates:
            r = _body_com_radius(b)
            if r > best_r:
                best_r = r
                ref_body = b
        if ref_body is None:
            raise RuntimeError(f"OB planes: MASTER {master_prefix} could not select reference body.")

        edges = [ref_body.edges.item(i) for i in range(ref_body.edges.count)]
        outer, inner, _ = _split_edges_by_radius(edges)
        chosen = outer if SIDE_CHOICE.lower() == "outer" else inner
        if not chosen:
            chosen = [(0.0, e) for e in edges]

        scored = []
        for _, e in chosen:
            L = _edge_length(e)
            Z = _edge_z_span(e)
            scored.append((LEN_WEIGHT * L + ZSPAN_WEIGHT * Z, e))
        scored.sort(key=lambda t: t[0], reverse=True)

        master_pts = master_cum = master_total = master_planes = None
        for _, e in scored[:min(EDGE_TRIES, len(scored))]:
            try:
                pts = _sample_edge_points(e, SPLINE_SAMPLE_POINTS)
                path = _make_root_path_from_3d_spline(root, pts)

                planes_for_master = []
                for f in sorted(FRACTIONS, reverse=True):
                    frac_key = int(round(f * 1000))
                    nm = OB_PLANE_NAME_FMT.format(prefix=f"OB{MASTER_OB_J}", frac=frac_key)
                    planes_for_master.append(_add_root_plane_along_path(root, path, f, nm))

                cum, total = _build_path_cumulative_lengths(pts)
                master_pts, master_cum, master_total, master_planes = pts, cum, total, planes_for_master
                break
            except Exception:
                continue

        if not master_planes:
            raise RuntimeError("OB planes: failed building MASTER planes.")

        path_data_by_j = {j: (master_pts, master_cum, master_total, master_planes) for (j, _) in obj_occs}
        step_timer.mark("10 Build OB planes (MASTER)")

        # 12.2) IB PLANES (MASTER=IB1) -> probably could be better if IB2
        ib_occ_root = None
        for i in range(blanket_final_occ2.childOccurrences.count):
            c = blanket_final_occ2.childOccurrences.item(i)
            if _safe_name(c.component).strip().lower() == "ib" or _safe_name(c).strip().lower() == "ib":
                ib_occ_root = c
                break
        if ib_occ_root is None:
            raise RuntimeError("Could not find 'IB' occurrence under BLANKET_FINAL.")

        ibj_occs = []
        for i in range(ib_occ_root.childOccurrences.count):
            o = ib_occ_root.childOccurrences.item(i)
            blob = (_safe_name(o) + " " + _safe_name(o.component)).lower()
            m = re.search(r"ib\s*([0-9]+)", blob)
            if m:
                ibj_occs.append((int(m.group(1)), o))
        ibj_occs.sort(key=lambda t: t[0])
        if not ibj_occs:
            raise RuntimeError("No IB{i} occurrences found under BLANKET_FINAL/IB.")

        _delete_existing_planes(root, r"^IB\d+_PLANE_\d{4}$")
        IB_FRACTIONS = sorted(FRACTIONS, reverse=True)

        MASTER_IB_J = 1
        ib_occ_master = None
        for (j, occj) in ibj_occs:
            if j == MASTER_IB_J:
                ib_occ_master = occj
                break
        if ib_occ_master is None:
            ib_occ_master = ibj_occs[0][1]
            MASTER_IB_J = ibj_occs[0][0]

        master_prefix = f"IB{MASTER_IB_J}"
        bodies_with_occ = _collect_solid_bodies_with_occ(ib_occ_master)
        nmax = _find_max_layer_index_for_prefix(master_prefix, bodies_with_occ)
        if nmax < 0:
            raise RuntimeError(f"IB planes: MASTER {master_prefix} has no Layer_N bodies.")
        ref_body, _ = _find_layer_body_with_occ(master_prefix, nmax, bodies_with_occ)
        if ref_body is None:
            raise RuntimeError(f"IB planes: MASTER {master_prefix} missing Layer_{nmax}.")

        e1, _ = _largest_edge_on_body(ref_body)
        e2, _ = _largest_connected_edge(ref_body, e1)
        if not e1 or not e2:
            raise RuntimeError("IB planes: could not find 2 connected edges.")

        pts = _sample_two_connected_edges_points(e1, e2, SPLINE_SAMPLE_POINTS)
        root_path = _make_root_path_from_3d_spline(root, pts)
        cum, total = _build_path_cumulative_lengths(pts)

        planes_for_master = []
        for f in IB_FRACTIONS:
            frac_key = int(round(f * 1000))
            nm = IB_PLANE_NAME_FMT.format(prefix=f"IB{MASTER_IB_J}", frac=frac_key)
            planes_for_master.append(_add_root_plane_along_path(root, root_path, f, nm))

        path_data_ib_by_j = {j: (pts, cum, total, planes_for_master) for (j, _) in ibj_occs}
        step_timer.mark("11 Build IB planes (MASTER)")

        # 13) Split OB/IB/VV by planes
        meta = {}
        TYPE_ORDER = {"armor": 0, "first_wall": 1, "layer": 2}

        # tag OB bodies
        for (j, ob_occ_j) in obj_occs:
            prefix = f"OB{j}".lower()
            for bb in list(ob_occ_j.bRepBodies):
                try:
                    if not bb.isSolid:
                        continue
                except Exception:
                    continue
                nm = _safe_name(bb).strip().lower()
                t = _native_token(bb)
                if not t:
                    continue
                if nm == f"{prefix}_armor":
                    meta[t] = (j, "armor", 0, "OB")
                elif nm == f"{prefix}_first_wall":
                    meta[t] = (j, "first_wall", 0, "OB")
                else:
                    m = re.match(rf"^{prefix}_layer_(\d+)$", nm)
                    if m:
                        meta[t] = (j, "layer", int(m.group(1)), "OB")

        # tag IB bodies
        for (j, ib_occ_j) in ibj_occs:
            prefix = f"IB{j}".lower()
            for (bb, _) in _collect_solid_bodies_with_occ(ib_occ_j):
                try:
                    if not bb.isSolid:
                        continue
                except Exception:
                    continue
                nm = _safe_name(bb).strip().lower()
                t = _native_token(bb)
                if not t:
                    continue
                if nm == f"{prefix}_armor":
                    meta[t] = (j, "armor", 0, "IB")
                elif nm == f"{prefix}_first_wall":
                    meta[t] = (j, "first_wall", 0, "IB")
                else:
                    m = re.match(rf"^{prefix}_layer_(\d+)$", nm)
                    if m:
                        meta[t] = (j, "layer", int(m.group(1)), "IB")

        # split OB
        split_ok = split_skip = 0
        for (j, ob_occ_j) in obj_occs:
            _, _, _, planes_for_j = path_data_by_j[j]
            planes_sorted = sorted([_as_native(pl) for pl in planes_for_j], key=_plane_frac_from_name, reverse=True)
            for pl_root_native in planes_sorted:
                try:
                    pl_for_ob = pl_root_native.createForAssemblyContext(ob_occ_j)
                except Exception:
                    pl_for_ob = pl_root_native

                bodies_now = [bb for bb in list(ob_occ_j.bRepBodies) if getattr(bb, "isSolid", True)]
                for bb in bodies_now:
                    t = _native_token(bb)
                    if not t or t not in meta:
                        continue
                    j0, kind0, n0, side0 = meta[t]
                    if side0 != "OB":
                        continue
                    try:
                        inp = split_feats_root.createInput(bb, pl_for_ob, True)
                        feat = split_feats_root.add(inp)

                        new_bodies = []
                        for ii in range(feat.bodies.count):
                            nb = feat.bodies.item(ii)
                            try:
                                if nb.isSolid:
                                    new_bodies.append(nb)
                            except Exception:
                                pass
                        if not new_bodies:
                            split_skip += 1
                            continue

                        del meta[t]
                        for nb in new_bodies:
                            nt = _native_token(nb)
                            if nt:
                                meta[nt] = (j0, kind0, n0, "OB")
                        split_ok += 1
                    except Exception:
                        split_skip += 1

        # split IB
        ib_split_ok = ib_split_skip = 0
        for (j, ib_occ_j) in ibj_occs:
            _, _, _, planes_for_j = path_data_ib_by_j[j]
            planes_sorted = sorted([_as_native(pl) for pl in planes_for_j], key=_plane_frac_from_name, reverse=True)

            for pl_root_native in planes_sorted:
                bodies_with_occ = _collect_solid_bodies_with_occ(ib_occ_j)
                for (bb, bb_occ) in bodies_with_occ:
                    t = _native_token(bb)
                    if not t or t not in meta:
                        continue
                    j0, kind0, n0, side0 = meta[t]
                    if side0 != "IB":
                        continue
                    try:
                        try:
                            pl_for_ib = pl_root_native.createForAssemblyContext(bb_occ)
                        except Exception:
                            pl_for_ib = pl_root_native

                        inp = split_feats_root.createInput(bb, pl_for_ib, True)
                        feat = split_feats_root.add(inp)

                        new_bodies = []
                        for ii in range(feat.bodies.count):
                            nb = feat.bodies.item(ii)
                            try:
                                if nb.isSolid:
                                    new_bodies.append(nb)
                            except Exception:
                                pass
                        if not new_bodies:
                            ib_split_skip += 1
                            continue

                        del meta[t]
                        for nb in new_bodies:
                            nt = _native_token(nb)
                            if nt:
                                meta[nt] = (j0, kind0, n0, "IB")
                        ib_split_ok += 1
                    except Exception:
                        ib_split_skip += 1

        # ======================================================
        # split VV by planes (robust: metadata-tracked)
        # ======================================================

        vv_plane_ok_ob = 0
        vv_plane_skip_ob = 0
        vv_plane_ok_ib = 0
        vv_plane_skip_ib = 0
        vv_plane_first_err = None

        vv_occ2 = _find_first_occurrence_by_hint(root, VV_HINT)
        if vv_occ2:

            def _vv_bodies_matching_prefix(prefix_lower: str):
                out = []
                for b in list(_solid_proxy_bodies(vv_occ2)):
                    nm = _safe_name(b).strip().lower()
                    if nm.startswith(prefix_lower):
                        out.append(b)
                return out

            def _split_vv_prefix_by_plane_list(prefix_lower: str, planes_native_sorted):
                nonlocal vv_plane_first_err
                ok = 0
                skip = 0

                for pl_root_native in planes_native_sorted:
                    bodies_now = _vv_bodies_matching_prefix(prefix_lower)
                    if not bodies_now:
                        continue

                    for vv_body in bodies_now:
                        try:
                            pl_for_vv = _plane_in_body_context(pl_root_native, vv_body)
                            inp = split_feats_root.createInput(vv_body, pl_for_vv, True)
                            split_feats_root.add(inp)
                            ok += 1
                        except Exception as e:
                            skip += 1
                            if vv_plane_first_err is None:
                                vv_plane_first_err = (
                                    f"VV split failed for {prefix_lower} with plane '{_safe_name(pl_root_native)}': {e}"
                                )

                return ok, skip

            # OB planes -> split VV_OB1..3 families
            for j in sorted(path_data_by_j.keys()):
                _, _, _, planes_for_j = path_data_by_j[j]
                planes_sorted = sorted([_as_native(pl) for pl in planes_for_j],
                                    key=_plane_frac_from_name, reverse=True)
                ok, sk = _split_vv_prefix_by_plane_list(f"vv_ob{j}", planes_sorted)
                vv_plane_ok_ob += ok
                vv_plane_skip_ob += sk

            # IB planes -> split VV_IB1..2 families
            for j in sorted(path_data_ib_by_j.keys()):
                _, _, _, planes_for_j = path_data_ib_by_j[j]
                planes_sorted = sorted([_as_native(pl) for pl in planes_for_j],
                                    key=_plane_frac_from_name, reverse=True)
                ok, sk = _split_vv_prefix_by_plane_list(f"vv_ib{j}", planes_sorted)
                vv_plane_ok_ib += ok
                vv_plane_skip_ib += sk

        step_timer.mark("13 Split OB/IB/VV by planes")

        # 14) FINALIZE: OB_SORTED / IB_SORTED + copy VV into buckets
        bf_comp = blanket_final_occ2.component

        def unique_child_name(parent_comp, base: str):
            existing = set()
            for i in range(parent_comp.occurrences.count):
                o = parent_comp.occurrences.item(i)
                try:
                    existing.add((_safe_name(o.component)).strip().lower())
                except Exception:
                    pass
            if base.lower() not in existing:
                return base
            idx = 2
            while f"{base}_{idx}".lower() in existing:
                idx += 1
            return f"{base}_{idx}"

        # OB_SORTED
        ob_sorted_name = unique_child_name(bf_comp, "OB_SORTED")
        _, ob_sorted_comp = new_component(bf_comp, ob_sorted_name)
        ob_j_comp = {}
        ob_k_comp = {}
        ob_items = []

        for (j, ob_occ_j) in obj_occs:
            pts, cum, total, _ = path_data_by_j[j]
            for bb in list(ob_occ_j.bRepBodies):
                try:
                    if not bb.isSolid:
                        continue
                except Exception:
                    continue
                t = _native_token(bb)
                if not t or t not in meta:
                    continue
                j0, kind0, n0, side0 = meta[t]
                if side0 != "OB":
                    continue
                try:
                    com = bb.physicalProperties.centerOfMass
                    k = _frac_to_k(_centroid_to_path_fraction(com, pts, cum, total))
                except Exception:
                    k = 1
                ob_items.append((j0, k, TYPE_ORDER.get(kind0, 9), n0, kind0, bb))

        ob_items.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
        ob_originals = []

        for (j0, k, _, n0, kind0, bb) in ob_items:
            if j0 not in ob_j_comp:
                _, cj = new_component(ob_sorted_comp, f"OB{j0}")
                ob_j_comp[j0] = cj
            key = (j0, k)
            if key not in ob_k_comp:
                _, ck = new_component(ob_j_comp[j0], f"k{str(k).zfill(2)}")
                ob_k_comp[key] = ck

            if kind0 == "armor":
                new_name = f"OB{j0}_{k}_armor"
            elif kind0 == "first_wall":
                new_name = f"OB{j0}_{k}_first_wall"
            else:
                new_name = f"OB{j0}_{k}_layer_{n0}"

            copy_body_into_def(ob_k_comp[key], bb, new_name)
            ob_originals.append(bb)

        # IB_SORTED
        ib_sorted_name = unique_child_name(bf_comp, "IB_SORTED")
        _, ib_sorted_comp = new_component(bf_comp, ib_sorted_name)
        ib_j_comp = {}
        ib_k_comp = {}
        ib_items = []

        for (j, ib_occ_j) in ibj_occs:
            pts, cum, total, _ = path_data_ib_by_j[j]
            for (bb, _) in _collect_solid_bodies_with_occ(ib_occ_j):
                try:
                    if not bb.isSolid:
                        continue
                except Exception:
                    continue
                t = _native_token(bb)
                if not t or t not in meta:
                    continue
                j0, kind0, n0, side0 = meta[t]
                if side0 != "IB":
                    continue
                try:
                    com = bb.physicalProperties.centerOfMass
                    k = _frac_to_k(_centroid_to_path_fraction(com, pts, cum, total))
                except Exception:
                    k = 1
                ib_items.append((j0, k, TYPE_ORDER.get(kind0, 9), n0, kind0, bb))

        ib_items.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
        ib_originals = []

        for (j0, k, _, n0, kind0, bb) in ib_items:
            if j0 not in ib_j_comp:
                _, cj = new_component(ib_sorted_comp, f"IB{j0}")
                ib_j_comp[j0] = cj
            key = (j0, k)
            if key not in ib_k_comp:
                _, ck = new_component(ib_j_comp[j0], f"k{str(k).zfill(2)}")
                ib_k_comp[key] = ck

            if kind0 == "armor":
                new_name = f"IB{j0}_{k}_armor"
            elif kind0 == "first_wall":
                new_name = f"IB{j0}_{k}_first_wall"
            else:
                new_name = f"IB{j0}_{k}_layer_{n0}"

            copy_body_into_def(ib_k_comp[key], bb, new_name)
            ib_originals.append(bb)

        if FINALIZE_DELETE_ORIGINALS:
            for bb in ob_originals + ib_originals:
                try:
                    (bb.nativeObject if bb.nativeObject else bb).deleteMe()
                except Exception:
                    pass
        elif FINALIZE_HIDE_ORIGINALS:
            for bb in ob_originals + ib_originals:
                try:
                    bb.isVisible = False
                except Exception:
                    pass

        # Copy VV into k buckets
        vv_copy_msg = "copy_vv_k: skipped (VV not found)."
        vv_occ3 = _find_first_occurrence_by_hint(root, VV_HINT)
        if vv_occ3:
            vv_copied, vv_hidden, vv_skipped = copy_vv_into_sorted(
                vv_occ=vv_occ3,
                ob_sorted_comp=ob_sorted_comp,
                ib_sorted_comp=ib_sorted_comp,
                path_data_by_j=path_data_by_j,
                path_data_ib_by_j=path_data_ib_by_j
            )
            vv_copy_msg = f"copy_vv_k: copied={vv_copied}, hidden={vv_hidden}, skipped={vv_skipped}."
        step_timer.mark("14 Finalize + copy VV")

        #15) 
        # 15.1) Port filling (fixed swap)
        if DO_PORT_FILLING and vv_occ3:
            pf_ok, pf_msg, _ = create_port_filling_block(
                root_comp=root,
                vv_occ=vv_occ3,
                port_body_name=PORT_FILL_PORTNAME,
                theta_deg=PORT_FILL_THETA_DEG,
                width_cm=PORT_FILL_WIDTH_CM,
                height_cm=PORT_FILL_HEIGHT_CM,
                swap_wh=PORT_FILL_SWAP_WH,
                use_square=PORT_FILL_USE_SQUARE,
                depth_scale=PORT_FILL_DEPTH_SCALE
            )
            vv_report.append("----- PORT FILLING -----")
            vv_report.append(pf_msg)

        # 15.2) de-overlap
        if DO_DEOVERLAP_15_2:
            deok, demsg, _ = boolean_deoverlap_blanket_final(bf_comp, min_keep_vol_cm3=1e-6, tools_chunk=40)
            vv_report.append("----- 15.2 DE-OVERLAP -----")
            vv_report.append(demsg)

        # 15.3) BF x VV cut
        if DO_BF_X_VV_CUT_15_3:
            bf_occ = blanket_final_occ2
            vv_occ_for_cut = _find_first_occurrence_by_hint(root, VV_HINT)
            bxok, bxmsg, _ = boolean_cut_blanket_final_by_vacuum_vessel_visible(
                root_comp=root,
                bf_occ=bf_occ,
                vv_occ=vv_occ_for_cut,
                tools_chunk=40,
                min_keep_vol_cm3=1e-6,
                vv_name_exclude_prefixes=("vv_port_",),
                vv_name_exclude_exact=("vv_divertor",)
            )
            vv_report.append("----- 15.3 BF x VV CUT -----")
            vv_report.append(bxmsg)

        # overlap verify (definition-space)
        if DO_OVERLAP_VERIFY:
            ov, msg = detect_overlaps_in_component(bf_comp, min_overlap_cm3=1e-6, only_visible=True)
            vv_report.append(msg)
            if ov:
                vv_report.append("Top overlaps:")
                for a, b, v in ov[:25]:
                    vv_report.append(f"  {a} x {b} => {v:.6g} cm^3")

        step_timer.mark("15 Deoverlap/CUT/verify")

        # 16) Export STEP + Tokamak_inputs.json
        export_name = "EUDEMO_complete_2_20"

        # Build STEP output *file* path
        step_file = os.path.join(export_folder, export_name + ".step")
        step_file = _ensure_parent_folder_exists(step_file)  # creates export_folder if missing

        # JSON path (same folder)
        inputs_json_path = os.path.join(export_folder, "Tokamak_inputs.json")
        inputs_json_path = _ensure_parent_folder_exists(inputs_json_path)

        # Export STEP
        exp = design.exportManager
        step_opts2 = exp.createSTEPExportOptions(step_file)
        if hasattr(step_opts2, "exportAsSolids"):
            step_opts2.exportAsSolids = True
        exp.execute(step_opts2)

        # Write JSON once
        write_user_inputs_json(inputs_json_path, root, step_export_path=step_file)

        step_timer.mark("16 Export STEP + JSON")

        # 17) Final report
        lines = []
        lines.append("FULL PIPELINE COMPLETED.\n")
        lines.append(f"STEP in : {STEP_PATH}")
        lines.append(f"STEP out: {step_file}")
        lines.append(f"JSON out: {inputs_json_path}\n")
        lines.append("----- WEDGE -----")
        lines.append(f"Wedge: 0° .. {ANGLE_MAX_DEG}° | Split0 ok/skip={ok0}/{sk0} | SplitA ok/skip={okA}/{skA}")
        lines.append(f"Removed: total={removed} (low={removed_low}, high={removed_high}) remove_fail={remove_fail}\n")
        lines.append("----- LAYERING -----")
        lines.extend(["  " + s for s in ob_summary])
        lines.extend(["  " + s for s in ib_summary])
        lines.append(vv_copy_msg + "\n")
        lines.append("----- VV REPORT -----")
        lines.extend(vv_report)
        lines.append("\n----- TIMINGS -----")
        lines.append(step_timer.format_table())

        final_msg = "\n".join(lines)
        _log_to_text_palette(app, final_msg + "\n")
        try:
            ui.messageBox("Tokamak_test_2 completed.\n\nSee Text Commands window for full report.")
        except Exception:
            pass

    except Exception:
        if ui:
            ui.messageBox("Failed:\n{}".format(traceback.format_exc()))
        else:
            raise