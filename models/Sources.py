import numpy as np
from itertools import count
from classes.Nodes import Nodes
# from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse
    
class IndepedentVoltageSources:
    def __init__(self, name, vp_node, vn_node, amp_ph_ph_rms, phase_deg, frequency_hz):
        self.name = name
        self.vp_node = vp_node
        self.vn_node = vn_node
        self.amp_ph_ph_rms = amp_ph_ph_rms
        self.phase_deg = phase_deg
        self.frequency_hz = frequency_hz
        self.DC = False
        # You are welcome to / may be required to add additional class variables   

    def get_nom_voltage(self):
        return self.amp_ph_ph_rms
    
    def get_current_voltage(self, t):
        if (self.DC):
            return self.get_nom_voltage()
        else:
            f = self.frequency_hz
            return math.sqrt(2/3) * self.amp_ph_ph_rms * math.cos(2 * math.pi * f * t + math.radians(self.phase_deg))
    
    # Some suggested functions to implement, 
    def assign_node_indexes(self,):
        pass
        
    def stamp_sparse(self,):
        pass

    def stamp_dense(self,):
        pass

    def stamp_open(self,):
        pass
    
    def __str__(self):
        return "V-{}-{}".format(self.amp_ph_ph_rms, id(self))
    
    def __repr__(self):
        return self.__str__()
    

class IndependentCurrentSource:
    # default object represents an independent current source
    def __init__(self, name, ip_node, in_node, amps):
        self.name = name
        self.ip_node = ip_node
        self.in_node = in_node
        self.amps = amps
        self.ecm_type = ""
        self.ecm_val = 0

    # Some suggested functions to implement, 
    def assign_node_indexes(self,):
        pass
        
    def stamp_sparse(self,):
        pass

    def stamp_dense(self,):
        pass

    def stamp_open(self,):
        pass
       
    def __str__(self): 
        return "I-{}-{}".format(self.amps, id(self))
    
    def __repr__(self):
        return self.__str__()    
        
    def stamp_sparse(self,):
        pass

    def stamp_dense(self,):
        pass

    def stamp_open(self,):
        pass
