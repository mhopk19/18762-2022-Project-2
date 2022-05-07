from __future__ import division
from models.Buses import Buses
import numpy as np


class Slack:

    def __init__(self,
                 Bus,
                 Vset,
                 ang,
                 Pinit,
                 Qinit):
        """Initialize slack bus in the power grid.

        Args:
            Bus (int): the bus number corresponding to the slack bus.
            Vset (float): the voltage setpoint that the slack bus must remain fixed at.
            ang (float): the slack bus voltage angle that it remains fixed at.
            Pinit (float): the initial active power that the slack bus is supplying
            Qinit (float): the initial reactive power that the slack bus is supplying
        """
        # You will need to implement the remainder of the __init__ function yourself.
        
        self.Bus = Bus
        self.Vset = Vset
        self.ang = ang
        self.Pinit = Pinit
        self.Qinit = Qinit
        
        self.Vr_set = Vset*np.cos(ang*np.pi/180)
        self.Vi_set = Vset*np.sin(ang*np.pi/180)
        
        self.Ir_init = (-self.Vr_set*self.Pinit - self.Vi_set*self.Qinit)/(Vset**2)
        self.Ii_init = (-self.Vi_set*self.Pinit + self.Vi_set*self.Qinit)/(Vset**2)
        

    def assign_nodes(self):
        """Assign the additional slack bus nodes for a slack bus.

        Returns:
            None
        """
        self.node_Vr = Buses.bus_map[self.Bus].node_Vr
        self.node_Vi = Buses.bus_map[self.Bus].node_Vi
        
        # these are the nodes for the current ammeter rows
        self.node_Ir = Buses._node_index.__next__()
        self.node_Ii = Buses._node_index.__next__()
    
    def assign_lambda_nodes(self):
        self.lambda_Vr = Buses.bus_map[self.Bus].lambda_Vr
        self.lambda_Vi = Buses.bus_map[self.Bus].lambda_Vi
        
        # these are the nodes for the current ammeter rows
        self.lambda_Ir = Buses._node_index.__next__()
        self.lambda_Ii = Buses._node_index.__next__()
    
    
    def stamp(self,Y,J, x_offset = 0, y_offset = 0):
        xx,yy = [x_offset, y_offset]
        i_node_r = self.node_Ir
        i_node_i = self.node_Ii
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        # independent voltage source stamping
        Y[i_node_r + xx][v_node_r + yy] += 1
        Y[i_node_i + xx][v_node_i + yy] += 1
        J[i_node_r + xx] += self.Vr_set
        J[i_node_i + xx] += self.Vi_set
        #
        Y[v_node_r + xx][i_node_r + yy] += 1
        Y[v_node_i + xx][i_node_i + yy] += 1
        
        return Y,J  
    