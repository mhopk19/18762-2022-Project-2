import numpy as np
from models.Buses import Buses

def initialize(size_Y, bus, generator, slack, flat_start = False, feasibility = False):
    
    v = np.zeros(size_Y, dtype=np.float)
    if flat_start:
        for ele in bus:
            v[ele.node_Vr] = 1
            v[ele.node_Vi] = 0
        for ele in generator:
            v[ele.node_Q] += (ele.Qmax+ele.Qmin)/2
        # leave slack currents initialized as 0?
    else:
        for ele in bus:
            v[ele.node_Vr] = ele.Vr_init
            v[ele.node_Vi] = ele.Vi_init
        for ele in generator:
            v[ele.node_Q] += -ele.Qinit
        for ele in slack:
            v[ele.node_Ir] = ele.Ir_init
            v[ele.node_Ii] = ele.Ii_init
    
    
    
    if (feasibility):
        lambda_current_epsilon = 0.001
        for b in bus:
            v[b.lambda_Vr] = lambda_current_epsilon
            v[b.lambda_Vr] = lambda_current_epsilon
    
    """
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
    """
    return v