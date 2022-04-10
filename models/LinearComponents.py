import numpy as np
from itertools import count
from classes.Nodes import Nodes
# from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse

class Resistors:
    def __init__(self, name, from_node, to_node, r):
        self.name = name
        self.from_node = from_node
        self.to_node = to_node
        self.r = r
        # You are welcome to / may be required to add additional class variables   

    # Some suggested functions to implement, 
    def assign_node_indexes(self, from_node, to_node):
        self.from_node = from_node
        self.to_node = to_node
        
    def stamp_sparse(self,):
        pass
    
class IndependentVoltageSource:
    def __init__(self, name, vp_node, vn_node, amp_ph_ph_rms, phase_deg, frequency_hz):
        self.name = name
        self.vp_node = vp_node
        self.vn_node = vn_node
        self.amp_ph_ph_rms = amp_ph_ph_rms
        self.phase_deg = phase_deg
        self.frequency_hz = frequency_hz
        # You are welcome to / may be required to add additional class variables   

    # Some suggested functions to implement, 
    def assign_node_indexes(self,vp_node, vn_node):
        self.vp_node = vp_node
        self.vn_node = vn_node        
        
    def stamp_sparse(self,):
        pass

    def stamp_dense(self,):
        pass

    def stamp_open(self,):
        pass