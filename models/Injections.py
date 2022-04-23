from __future__ import division
from itertools import count
from models.Buses import Buses
import math
import scripts.global_vars as gv

class Injections:
    _ids = count(0)
    base = gv.global_vars.MVAbase

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
        
    
    def assign_nodes(self):
        """Assign the additional bus nodes for injections

        Returns:
            None
        """
        # these are the nodes for the current ammeter rows
        self.node_R = Buses._node_index.__next__()
        self.node_I = Buses._node_index.__next__()
        
    
    def stamp(self, Y, J, prev_v):
            
        return Y, J
        
        
        
        
    