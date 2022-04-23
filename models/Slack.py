from __future__ import division
from models.Buses import Buses


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

        # initialize nodes
        self.node_Vr_Slack = None
        self.node_Vi_Slack = None
        
        self.Bus = Bus
        self.Vset = Vset
        self.ang = ang
        self.Pinit = Pinit
        self.Qinit = Qinit

    def assign_nodes(self):
        """Assign the additional slack bus nodes for a slack bus.

        Returns:
            None
        """
        # these are the nodes for the current ammeter rows
        self.node_Vr_Slack = Buses._node_index.__next__()
        self.node_Vi_Slack = Buses._node_index.__next__()
    
    
    def stamp(self,Y,J, x_offset = 0, y_offset = 0):
        xx,yy = [x_offset, y_offset]
        i_node_r = self.node_Vr_Slack
        i_node_i = self.node_Vi_Slack
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        # independent voltage source stamping
        Y[i_node_r + xx][v_node_r + yy] = 1
        Y[i_node_i + xx][v_node_i + yy] = 1
        J[i_node_r + xx] = self.Vset
        J[i_node_i + xx] = 0
        #
        Y[v_node_r + xx][i_node_r + yy] = 1
        Y[v_node_i + xx][i_node_i + yy] = 1
        J[v_node_r + xx] = 0
        J[v_node_i + xx] = 0
        
        return Y,J  
    