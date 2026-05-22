# ===========================================
# TEA: Manufacturing Cost Model
# per equations 3.1–3.9 in "Costing Designs"
# ===========================================

from typing import List, Dict
from cost_variables import material_cost_database

import sys
import os
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'materials'))
sys.path.append(module_path)
import materials

import openmc

# --------------------------
# Material Class
# --------------------------
class ManufacturingMaterial:
    """
    Represents a material with base unit cost (Cmt)
    and process-suitability coefficients (Cmp).
    """

    def __init__(self, material: openmc.Material, price_override: float, cmp_map: Dict[str, float]):
        """
        Parameters
        ----------
        material : openmc.Material class
            Nuclide types can be accessed via material.nuclides[i].name. 
            Nuclide fraction in the material can be accessed via material.nuclides[i].percent. 
            Nuclide percent type can be accessed via material.nuclides[i].percent_type. 
        price_override : float
            Cost per mm³ (USD/mm³), optional override for computed cost
        cmp_map : dict
            Dictionary mapping process names to relative suitability (Cmp)

        TODO: 
            - add imppurity control cost (e.g. super conuctor magnets?), toxic material handling cost
            - allow price overide to be in $/kg as well as $/mm³
        """
        self.material = material
        self.name = material.name
        self.density = material.density                                                                    # g/cm³
        self.cmp_map = cmp_map
        self.price_override = price_override
        self.price_volume, self.price_mass = self._compute_price()                                         # USD/mm³, USD/kg
        

    def _compute_price(self) -> float:
        """
        Get base material cost (Cmt).
        """
        print(f"Calculating material cost for {self.name}...")
        if self.price_override is not None:
            cost_v = self.price_override                                                                   # USD/mm³
            cost_m = cost_v * 1000 * 1000 / self.density                                                   # USD/kg---[USD/mm³] [cm³ / g] * [1000 mm³ / cm³] [1000 g / kg]
            print(f"\tMaterial cost (Cmt) from override: ${cost_m:.2f}/kg")
            print(f"\tMaterial cost (Cmt) from override: ${cost_v:.6f}/mm³")
            return cost_v, cost_m
        else:
            material_masses = {}
            for nuclide in self.material.nuclides:
                # Extract atomic symbol (letters only) from nuclide name
                symbol = ''.join(filter(str.isalpha, nuclide.name))
                if nuclide.percent_type == 'wo':
                    mass = nuclide.percent/100
                else:
                    raise ValueError(f"Unsupported percent_type '{nuclide.percent_type}' for nuclide {nuclide.name}. Only 'wo' is supported.")
                material_masses[symbol] = material_masses.get(symbol, 0.0) + mass

            # Calculate total cost
            cost = 0.0                                                                                     # USD/kg
            for symbol, mass in material_masses.items():
                if symbol in material_cost_database:
                    unit_cost = material_cost_database[symbol].cost * mass
                    cost += unit_cost
                else:
                    print(f"\tWarning: Component, {symbol}, not found in material cost database.")
            print(f"\tMaterial cost (Cmt): ${cost:.2f}/kg")
            cost_m = cost                                                                                  # USD/kg                                    
            # convert USD/kg to USD/mm³ using density [g/cm³]
            # 1 mm³ = 1e-3 cm³ and 1 kg = 1000 g, so
            # cost_v = cost_m * density [g/cm³] * (1 kg / 1000 g) * (1 cm³ / 1000 mm³)
            cost_v = cost * self.density / 1000 / 1000                                                     # USD/mm³
            print(f"\tMaterial cost (Cmt): ${cost_v:.6f}/mm³")
            return cost_v, cost_m


    def get_cmp(self, process_name: str) -> float:
        """
        Get material–process suitability coefficient (Cmp).
        Default = 1.5 if process not listed.

        TODO: Add a warning that process is not listed and a defult is being used
        """
        return self.cmp_map.get(process_name, 1.5)

    def __repr__(self):
        return f"<Material: {self.name}>"


# --------------------------
# ManufacturingProcess Class
# --------------------------
class ManufacturingProcess:
    """
    Represents a manufacturing process with cost and performance data.
    """

    def __init__(self, name: str, alphaT: float, beta: float, description: str = ""):
        self.name = name
        self.alphaT = alphaT                                                                            # Product of operating and overhead cost per sec (α) and time (T)
        self.beta = beta                                                                                # Tooling cost for ideal design (β)
        self.description = description

    def basic_processing_cost(self, production_qty: int) -> float:
        """Equation [3.3]: Pc = α*T + β/N"""
        print(f"\tCalculating basic processing cost for {self.name}...")
        Pc = self.alphaT + self.beta / production_qty
        print(f"\t\tIdeal processing cost (Pc): ${Pc:.2f} per component")
        return Pc

    def __repr__(self):
        return f"<Process: {self.name}>"


# --------------------------
# Component Class
# --------------------------
class ManufacturingComponent:
    """
    Represents a component that can have multiple sequential processes.
    """

    def __init__(self, name: str, volume_mm3: float, material: ManufacturingMaterial,
                 processes: List[ManufacturingProcess], production_qty: int,
                 cc=1.0, cs=1.0, ct=1.0, cf=1.0, waste_coeff=1.0):
        """
        Parameters
        ----------
        name : str
            Name of the component
        volume_mm3 : float
            Volume of the component in mm³
        material : Material
            Material object
        processes : list of ManufacturingProcess objects
            List of manufacturing processes applied to the component
        production_qty : int
            Production quantity (N)
        cc : float
            Complexity coefficient (Cc)
        cs : float
            Size coefficient (Cs)
        ct : float
            Tolerance coefficient (Ct)
        cf : float
            Finish coefficient (Cf)
        waste_coeff : float
            Waste coefficient (Wc)
        
        
        """
        self.name = name
        self.Vf = volume_mm3
        self.material = material
        self.processes = processes
        self.N = production_qty
        self.Cc = cc
        self.Cs = cs
        self.Ct = ct
        self.Cf = cf
        self.Wc = waste_coeff

    # --- Cost Model Equations ---

    def material_cost(self) -> float:
        """Eq. [3.6]–[3.7]: Mc = V * Cmt"""
        V_input = self.Vf * self.Wc                                                             # mm³---[mm³] [ ]
        Mc = V_input * self.material.price_volume                                               # USD---[mm³] [USD/mm³]
        print(f"\tCalculating material cost for component {self.name}...")
        print(f"\t\tMaterial cost (MC): {Mc:.2f} per component")
        return Mc

    def relative_cost_coefficient(self, process: ManufacturingProcess) -> float:
        """
        Parameters
        ----------
        process :  ManufacturingProcess object
            Name of a process (e.g., CNC)
        
        Returns
        -------
        Rc : float
            Relative cost coefficient---Rc = Cmp * Cc * Cs * Cft (Eq. 3.5) 
        """
        print(f"\tCalculating relative cost coefficient for processing {self.material.name} by {process.name} to make {self.name}...")
        cmp = self.material.get_cmp(process.name) 
        Cft = max(self.Ct, self.Cf)
        Rc = cmp * self.Cc * self.Cs * Cft
        print(f"\t\tRelative cost coefficient (Rc): {Rc:.2f}")
        return Rc

    def manufacturing_cost(self) -> float:
        """Eq. [3.2]: Mi = Mc + Σ(Rc_i * Pc_i)"""
        print(f"Calculating manufacturing cost for component {self.name}...")
        Mc = self.material_cost()
        total_process_cost = 0
        for p in self.processes:
            Pc = p.basic_processing_cost(self.N)
            Rc = self.relative_cost_coefficient(p)
            total_process_cost += Rc * Pc
            total_cost = Mc + total_process_cost
            mass = self.material.density * self.Vf / 1000 / 1000     # kg---[g/cm³] [mm³] * [1 kg/1000 g] [1 cm³/1000 mm³]
            total_unit_cost = total_cost / mass                      # $/kg
        print(f"\tTotal manufacturing cost (Mi): ${total_cost:.2f} per component")
        print(f"\tTotal manufacturing cost (Mi): ${total_unit_cost:.2f} per kg")
        return total_unit_cost

    def summary(self):
        summary_data = {
            "Component": self.name,
            "Material": self.material.name,
            "Material cost (Mc)": self.material_cost(),
            "Processes": [],
            "Total manufacturing cost (Mi)": self.manufacturing_cost()
        }
        for p in self.processes:
            Pc = p.basic_processing_cost(self.N)
            Rc = self.relative_cost_coefficient(p)
            summary_data["Processes"].append({
                "Process": p.name,
                "Pc": Pc,
                "Rc": Rc,
                "Rc*Pc": Rc * Pc
            })
        return summary_data


# -------------------------------------------
# Build Processes (Fig. )
# -------------------------------------------
def build_processes():
    data = {
        "CMC":  (22.29, 5886.684),
        "IC":   (16.74, 5942.739),
        "SC":   (12.87, 5981.826),
        "SM":   (9.90, 6011.823),
        "CPM":  (1.80, 5220.0576),
        "GDC":  (1.19, 79778.7694),
        "IM":   (0.80, 99992.0),
        "PDC":  (0.80, 99992.0),
        "VF":   (11.58, 2231.3205),
        "PM":   (1.98, 55938.0458),
        "CDF":  (0.93, 20640.3311),
        "CF":   (0.55, 25953.1792),
        "CH":   (0.50, 17235.4812),
        "SMW":  (0.39, 12508.4652),
        "CNC":  (7.84, 2480.3779),
        "MM":   (5.28, 4945.5924),
        "AM":   (1.24, 50283.1242),
        "HCEM": (3.58, 29996.1574),
        "CCEM": (3.00, 24983.0664),
        "CEP":  (2.60, 99503.685),
        "CM5":  (1927.33, 2451.4918),
        "CM2.5":(957.46, 2567.9385)
    }

    return {name: ManufacturingProcess(name=name, alphaT=a, beta=b)
            for name, (a, b) in data.items()}


# ------------------------------------------
# Build Materials (Fig. 3.22/3.7)
# ------------------------------------------

def build_materials():
    """
    Generic map for alloy steel:
    ----------------------------
    cmp_map={
        "AM": 2.5,
        "CCEM": 2,
        "CDF": 2,
        "CF": 2,
        "CH": 2,     
        "CM2.5": 1,   
        "CM5": 1,
        "CMC": 1.3,
        "CNC": 2.5,
        "HCEM": 2,
        "IC": 1,
        "MM": 2.5,
        "PM": 1.1,
        "SM": 1.3,
        "SC": 1.3,
        "SMW": 1.5
    }
    Generic map for stainless steel:
    --------------------------------
    cmp_map={ 
        "AM": 4,
        "CCEM": 2,
        "CDF": 2,
        "CF": 2,
        "CH": 2,     
        "CM2.5": 1,   
        "CM5": 1,
        "CMC": 1.5,
        "CNC": 4,
        "HCEM": 2,
        "IC": 1,
        "MM": 4,
        "PM": 1.1,
        "SM": 1.5,
        "SC": 1.5,
        "SMW": 1.5 
    }
    Generic map for low carbon steel:
    ---------------------------------
    cmp_map={
        "AM": 1.4,
        "CCEM": 1.3,
        "CDF": 1,
        "CF": 1.3,
        "CH": 1.3,     
        "CM2.5": 1,   
        "CM5": 1,
        "CMC": 1.2,
        "CNC": 1.4,
        "HCEM": 1.3,
        "IC": 1,
        "MM": 1.4,
        "PM": 1.2,
        "SM": 1.2,
        "SC": 1.2,
        "SMW": 1.2    
    }
    Generic map for cast iron:
    --------------------------
    cmp_map={
        "AM": 1.2,     
        "CM2.5": 1,   
        "CM5": 1,
        "CMC": 1,
        "CNC": 1.2,
        "IC": 1,
        "MM": 1.2,
        "PM": 1.6,
        "SM": 1,
        "SC": 1    
    }
    """
    print('Adding Eurofer97 material...')
    eurofer97 = materials.eurofer97(8.0)   # density in g/cc
    eurofer97.name = "Eurofer97"
    eurofer97_MM = ManufacturingMaterial(
        material = eurofer97, 
        price_override = None,
        cmp_map={                          # based on generic stainless steel
            "AM": 4,
            "CCEM": 2,
            "CDF": 2,
            "CF": 2,
            "CH": 2,     
            "CM2.5": 1,   
            "CM5": 1,
            "CMC": 1.5,
            "CNC": 4,
            "HCEM": 2,
            "IC": 1,
            "MM": 4,
            "PM": 1.1,
            "SM": 1.5,
            "SC": 1.5,
            "SMW": 1.5    
        }
    )
    return {eurofer97_MM.name: eurofer97_MM}

ManufacturingMaterials_Database = build_materials()
ManufacturingProcesses_Database = build_processes()