#!/usr/bin/env python3

#import pyne
from pyne.material import Material as PyneMaterial
from openmc import Material, Materials
from xml.etree import ElementTree
import pyne.nucname as nucname

def make_openmc_mat(pyne_mat):
    xml_str = pyne_mat.openmc()
    xml_elem = ElementTree.fromstring(xml_str)
    openmc_mat = Material.from_xml_element(xml_elem)
    return openmc_mat

def make_materials_xml():

    steel_composition = {
        'Fe' : 98.861,
        'C' : 0.360,
        'S' : 0.042,
        'Si': 0.061,
        'P' : 0.056,
        'Mn': 0.620,
    }
    #carbon_nucids = [ nucname.id(nuclide) for nuclide in ['C12','C13'] ]
    steel = PyneMaterial(steel_composition)
    # Density in g/cc
    steel.density = 7.8
    steel = steel.expand_elements()
    steel.metadata["mat_number"]=1
    steel.metadata["name"]="SS"

    #Lithium lead
    lithium_lead_composition = {
        'Pb' : 0.83,
        'Li' : 0.17,
    }
    lithium_lead = PyneMaterial(lithium_lead_composition)
    # Density in g/cc
    lithium_lead.density = 1.593 # g/cc
    lithium_lead = lithium_lead.expand_elements()
    lithium_lead.metadata["mat_number"]=2
    lithium_lead.metadata["name"]="PbLi"

    tungsten_composition = { 'W' : 1.0 }
    tungsten = PyneMaterial(tungsten_composition)
    tungsten.density = 19.30 # g/cc
    tungsten.metadata["mat_number"]=3
    tungsten.metadata["name"]="W"

    pyne_mats = [ steel, lithium_lead, tungsten ]
    mat_list = [ make_openmc_mat(mat) for  mat in pyne_mats ]
    mats = Materials(mat_list)
    mats.export_to_xml()
    print("Created materials.xml")

    colors = ["green", "blue", "orange"]
    mat_dict = dict(zip(mat_list, colors))
    return mat_dict
