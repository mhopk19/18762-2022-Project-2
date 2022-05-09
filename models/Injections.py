from __future__ import division
from itertools import count
from models.Buses import Buses
import math
import scripts.global_vars as gv

class Injections:
    _ids = count(0)
    base = gv.global_vars.MVAbase
    init = False 
    node_neutral = None
    
    def __init__(self,
                 Bus):
        """Initialize an instance of a PQ or ZIP load in the power grid.

        Args:
            Bus (int): the bus where the load is located
            P (float): the active power of a constant power (PQ) load. [entered in MW]
            Q (float): the reactive power of a constant power (PQ) load. [entered in Mvar]
            IP (float): the active power component of a constant current load.
            IQ (float): the reactive power component of a constant current load.
            ZP (float): the active power component of a constant admittance load.
            ZQ (float): the reactive power component of a constant admittance load.
            area (int): location where the load is assigned to.
            status (bool): indicates if the load is in-service or out-of-service.
        """

        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.
        
        self.Bus = Bus
        
        # create the neutral node
        if (Injections.init == False):
            Injections.neutral_v = Buses._node_index.__next__()
            Injections.neutral_i = Buses._node_index.__next__()
            Injections.init = True
    
    def assign_nodes(self):
        """Assign the additional bus nodes for injections
    
        Returns:
            None
        """
        # these are the nodes for the current ammeter rows
        self.node_Vsr = Buses._node_index.__next__()
        self.node_Isr = Buses._node_index.__next__()
        self.node_Vsi = Buses._node_index.__next__()
        self.node_Isi = Buses._node_index.__next__()
        self.node_Vr = Buses._node_index.__next__()
        self.node_Ir = Buses._node_index.__next__()
        self.node_Vi = Buses._node_index.__next__()
        self.node_Ii = Buses._node_index.__next__()
        self.neutral_v = Injections.neutral_v
        self.neutral_i = Injections.neutral_i
     
    def stamp(self, Y, J, prev_v):
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        
        #J[v_node_r] = prev_lambda[self.Bus]
        #J[v_node_i] = prev_lambda[self.Bus + 1]
        
        # default resistance of 1
        
        # grounded source nodes
        Y[self.node_Vsr][self.node_Isr] = 1
        Y[self.node_Isr][self.node_Vsr] = 1

        Y[self.node_Vsi][self.node_Isi] = 1
        Y[self.node_Isi][self.node_Vsi] = 1
        
        
        # neutral node
        Y[self.neutral_v][self.neutral_i] = 1
        Y[self.neutral_i][self.neutral_v] = 1
        
        # independent real feasibility voltage
        Y[self.node_Vr][self.node_Ir] = 1
        Y[self.node_Ir][self.node_Vr] = 1
        
        J[self.node_Ir] = 0#prev_v[self.node_Vr]
        
        # lambda = current identity
        Y[self.node_Ir][Buses.bus_map[self.Bus].lambda_Vr] = -1
        
        # independent imag feasibility voltage
        Y[self.node_Vi][self.node_Ii] = 1
        Y[self.node_Ii][self.node_Vi] = 1
        
        J[self.node_Ii] = 0#prev_v[self.node_Vi]
        
        # lambda = current identity
        Y[self.node_Ii][Buses.bus_map[self.Bus].lambda_Vi] = -1
        
        
        # voltage controlled current source for real circuit
        Y[self.node_Vsr][self.node_Vr] = 1
        Y[self.node_Vsr][self.neutral_v] = -1
        Y[v_node_r][self.neutral_v] = -1
        Y[v_node_r][self.node_Vr] = 1
        
        # voltage controlled current source for imag circuit
        Y[self.node_Vsi][self.node_Vi] = 1
        Y[self.node_Vsi][self.neutral_v] = -1
        Y[v_node_i][self.neutral_v] = -1
        Y[v_node_i][self.node_Vi] = 1
        
        
        
        return Y, J
        
        
        
        
    