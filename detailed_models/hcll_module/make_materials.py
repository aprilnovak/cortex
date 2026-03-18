#!/usr/bin/env python3

from openmc import Material, Materials
from xml.etree import ElementTree
import materials as fusion_mats

def make_materials_xml():
    steel = fusion_mats.ss316(7.8)
    steel.id=1
    steel.name="SS"

    #Lithium lead
    lithium_lead = fusion_mats.PbLi(0.6, 1.593)
    lithium_lead.id=2
    lithium_lead.name="PbLi"

    tungsten = fusion_mats.W(19.30)
    tungsten.id=3
    tungsten.name="W"

    pyne_mats = [ steel, lithium_lead, tungsten ]
    mats = Materials(pyne_mats)
    mats.export_to_xml()
    print("Created materials.xml")

    colors = ["green", "blue", "orange"]
    mat_dict = dict(zip(pyne_mats, colors))
    return mat_dict
