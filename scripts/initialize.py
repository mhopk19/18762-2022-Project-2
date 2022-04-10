import numpy as np
from models.Buses import Buses

def initialize(size_Y, voltage_components):
    # initial conditions for GS-4
    #v = [1., 0., 1., 0., 1., 0., 1., 0., 0., 0., 0.]
    
    v = np.zeros(size_Y)
    
    # FLAT START conditions
    for v_comp in voltage_components:
        v_node_r = Buses.bus_map[v_comp.Bus].node_Vr
        v_node_i = Buses.bus_map[v_comp.Bus].node_Vi
        v[v_node_r] = 1
        v[v_node_i] = 0
    
    return v