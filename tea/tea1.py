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

        n = len(processes)

        self.Cc = Cc if Cc is not None else [1.0] * n
        self.Cs = Cs if Cs is not None else [1.0] * n
        self.Ct = Ct if Ct is not None else [1.0] * n
        self.Cf = Cf if Cf is not None else [1.0] * n
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
            Cmp = self.material.get_Cmp(p.name)

            Rc = Cmp * self.Cc[i] * self.Cs[i] * self.Ct[i] * self.Cf[i]

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