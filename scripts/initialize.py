import numpy as np
from models.Buses import Buses

def initialize(size_Y, voltage_components, bus, feasibility = False):
    # initial conditions for GS-4
    #v = [1., 0., 1., 0., 1., 0., 1., 0., 0., 0., 0.]  
    
    v = np.zeros(size_Y)
    
    # FLAT START conditions
    for v_comp in voltage_components:
        v_node_r = Buses.bus_map[v_comp.Bus].node_Vr
        v_node_i = Buses.bus_map[v_comp.Bus].node_Vi
        v[v_node_r] = 1
        v[v_node_i] = 0
        
    slack_current_epsilon = 0.00001
    
    if (feasibility):
        n_bus = 0
        for b in bus:
            n_bus = n_bus + 1
        
        n_comps = size_Y - (2 * n_bus)
        
        for i in range(n_comps,size_Y):
            v[i] = slack_current_epsilon
    
    return v