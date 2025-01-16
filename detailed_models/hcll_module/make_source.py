#!/usr/bin/env python3
import openmc
from openmc.stats import (
    CartesianIndependent,
    Monodirectional,
    Discrete,
    Uniform,
    Point
)

def make_plane_source():
    # Specify spatial distribution for source location
    # Specify each axis direction independently
    x_disbn = Uniform(-23.1,23.1)
    y_disbn = Uniform(-83.5,83.5)
    z_disbn = Discrete(x=[69.5],p=[1.0])
    spatial_disbn = CartesianIndependent(x_disbn,y_disbn,z_disbn)

    # Specify an angular distribution for flux directionality
    angle_disbn = Monodirectional([0.,0.,-1.])

    # Specify an energy distribution: 14 MeV neutrons
    energy_disbn = Discrete(x=[14.0e6],p=[1.0])

    source = openmc.Source(space=spatial_disbn,
                           angle=angle_disbn,
                           energy=energy_disbn,
                           particle='neutron')
    return source
