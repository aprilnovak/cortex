"""
TEA and DFM equations.


########################################################################
DESCRIPTION
########################################################################
The 'Component' class contains the inputs that define a component (e.g. 
volume, material, fabrcation process, etc.). And the TEA and DFM 
functions that assess the cost to manufacture a component as defined.

########################################################################
CLASSES
########################################################################
------------------------------------------------------------------------
Component (name: str, volume_mm3: float, material: Material, ...)
    Represents a component. 
"""


from process_def import Process
from material_def import Material
from component_def import Geometry
from typing import List

class Component:
    """
    Represents a component. 
    
    
    ########################################################################
    DESCRIPTION
    ########################################################################
    The 'Component' class contains the inputs that define a component (e.g. 
    volume, material, fabrcation process, etc.). And the TEA and DFM 
    functions that assess the cost to manufacture a component as defined.

    ########################################################################
    ATTRIBUTES
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the component

    ------------------------------------------------------------------------
    volume_mm3 : float
        Volume of the component in mm³

    ------------------------------------------------------------------------
    material : Material
        'Material' object

    ------------------------------------------------------------------------
    processes : list of 'Process' objects
        List of processes applied to the component

    ------------------------------------------------------------------------
    production_qty : int
        Production quantity (N)

    ------------------------------------------------------------------------
    geometry : Geometry, optional
        'Geometry' object. When given, Cc, Cs, Ct, and Cf are looked up
        per process from the geometry (geometry.get_Cc(process.name),
        etc.) for any of Cc/Cs/Ct/Cf not explicitly provided below.

    ------------------------------------------------------------------------
    cc : float
        Complexity coefficient (Cc)

    ------------------------------------------------------------------------
    cs : float
        Size coefficient (Cs)

    ------------------------------------------------------------------------
    ct : float
        Tolerance coefficient (Ct)

    ------------------------------------------------------------------------
    cf : float
        Finish coefficient (Cf)

    ------------------------------------------------------------------------
    waste_coeff : float
        Waste coefficient (Wc)
    
    ########################################################################
    METHODS
    ########################################################################
    ------------------------------------------------------------------------
    material_cost(verbose=False)
        Calculate material cost from composition

    ------------------------------------------------------------------------
    get_mass_kg()
        Convert volume (mm³) to mass (kg) using density (g/cm³)

    ------------------------------------------------------------------------
    manufacturing_cost()
        Calculate total manufacturing cost for the component

    ------------------------------------------------------------------------
    summary()
        Generate a summary of the cost breakdown for the component
    """

    def __init__(self,
                 name: str,
                 volume_mm3: float,
                 material: Material,
                 processes: List[Process],
                 production_qty: int,
                 geometry: Geometry = None,
                 Cc = None,
                 Cs = None,
                 Ct = None,
                 Cf = None,
                 Wc = None):

        self.name = name
        self.Vf = volume_mm3
        self.material = material
        self.processes = processes
        self.N = production_qty
        self.geometry = geometry

        n = len(processes)

        if Cc is not None:
            self.Cc = Cc
        elif geometry is not None:
            self.Cc = [geometry.get_Cc(p.name) for p in processes]
        else:
            self.Cc = [1.0] * n

        if Cs is not None:
            self.Cs = Cs
        elif geometry is not None:
            self.Cs = [geometry.get_Cs(p.name) for p in processes]
        else:
            self.Cs = [1.0] * n

        if Ct is not None:
            self.Ct = Ct
        elif geometry is not None:
            self.Ct = [geometry.get_Ct(p.name) for p in processes]
        else:
            self.Ct = [1.0] * n

        if Cf is not None:
            self.Cf = Cf
        elif geometry is not None:
            self.Cf = [geometry.get_Cf(p.name) for p in processes]
        else:
            self.Cf = [1.0] * n

        self.Wc = Wc if Wc is not None else [1.0] * n

    def material_cost(self, verbose=False) -> float:
        """
        Calculate material cost from composition
        """

        mass_kg = self.get_mass_kg()              # kg

        total_Wc = 1.0
        for wc in self.Wc:
            total_Wc *= wc

        effective_mass = mass_kg * total_Wc       # include waste

        Mc = effective_mass * self.material.cost  # $

        if verbose:
            print(f"[{self.name}] material cost calculation:")
            print(f"  Component mass (kg): {mass_kg:.4f}")
            print(f"  Waste coefficient: {self.Wc}")
            print(f"  Material mass (kg): {effective_mass:.4f}")
            print(f"  Material cost ($/kg): {self.material.cost:.2f}")
            print(f"  Material cost ($): {Mc:.2f}")

        return Mc

    
    def get_mass_kg(self) -> float:
        """
        Convert volume (mm³) to mass (kg) using density (g/cm³)
        """

        return self.material.density * self.Vf / 1e6 # kg---[g/cm³] [mm³] * [1 kg/1000 g] [1 cm³/1000 mm³]

    def get_Rc(self, i: int) -> float:
        """
        Relative cost coefficient (Rc) for the i-th process:

        Rc = Cmp * Cc * Cs * max(Ct, Cf)

        Only the greater of the tolerance (Ct) and finish (Cf) coefficients
        is applied, since achieving the tighter of the two typically also
        satisfies the other.
        """

        p = self.processes[i]
        Cmp = self.material.get_Cmp(p.name)

        return Cmp * self.Cc[i] * self.Cs[i] * max(self.Ct[i], self.Cf[i])

    def fabrication_cost_curve(self, quantities):
        """
        Fabrication cost (Rc * Pc) for every process, evaluated at the given
        production quantity/quantities - both at this component's actual Rc
        and at the DFM-ideal Rc (=1.0, i.e. Cmp=Cc=Cs=Ct=Cf=1). Used for the
        "cost vs. production quantity" sensitivity view; quantities may be a
        scalar or an array (Process.basic_processing_cost broadcasts).


        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        quantities : int, float, or array-like
            Production quantity/quantities (N) to evaluate at

        ########################################################################
        RETURNS
        ########################################################################
        ------------------------------------------------------------------------
        rows : list of dict
            One entry per process: {'process', 'Rc_selected', 'cost_ideal', 'cost_selected'}
        """

        rows = []
        for i, p in enumerate(self.processes):
            Rc_selected = self.get_Rc(i)
            Pc = p.basic_processing_cost(quantities)
            rows.append({
                "process": p,
                "Rc_selected": Rc_selected,
                "cost_ideal": Pc,
                "cost_selected": Rc_selected * Pc,
            })
        return rows

    def manufacturing_cost(self, verbose=False) -> dict:
        """
        M = Mc + Σ(Rc_i * Pc_i)
        """

        if verbose:
            print(f"Calculating manufacturing cost for component {self.name}...")

        Mc = self.material_cost()

        process_costs = []
        breakdown = {"Material": Mc}
        total_process_cost = 0

        for i, p in enumerate(self.processes):

            Pc = p.basic_processing_cost(self.N)
            Rc = self.get_Rc(i)

            process_cost = Rc * Pc
            total_process_cost += process_cost
            
            process_costs.append({
                "Process": p.name,
                "Pc": Pc,
                "Rc": Rc,
                "Cost": process_cost
            })

            breakdown[p.name] = process_cost

            if verbose:
                print(f"\tProcess: {p.name}")
                print(f"\t  Pc: {Pc:.2f}")
                print(f"\t  Rc: {Rc:.3f}")
                print(f"\t  Cost: {process_cost:.2f}")

        total_cost = Mc + total_process_cost
        mass = self.get_mass_kg()
        total_unit_cost = total_cost / mass 

        total_Wc = 1.0
        for wc in self.Wc:
            total_Wc *= wc

        if verbose:
            print(f"\tTotal manufacturing cost: ${total_cost:.2f} per component")
            print(f"\tTotal manufacturing cost: ${total_unit_cost:.2f} per kg")

        return {
            "Component": self.name,
            "Mass": mass,
            "Wc": total_Wc,
            "Material": self.material.name,
            "Processes": process_costs,
            "Cost Breakdown": breakdown,
            "Processing cost": total_process_cost,
            "Material cost": Mc,
            "Total cost": total_cost,
            "Unit cost": total_unit_cost
        }