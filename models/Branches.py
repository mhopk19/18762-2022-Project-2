from __future__ import division
from itertools import count
from models.Buses import Buses

class Branches:
    _ids = count(0)

    def __init__(self,
                 from_bus,
                 to_bus,
                 r,
                 x,
                 b,
                 status,
                 rateA,
                 rateB,
                 rateC):
        """Initialize a branch in the power grid.

        Args:
            from_bus (int): the bus number at the sending end of the branch.
            to_bus (int): the bus number at the receiving end of the branch.
            r (float): the branch resistance
            x (float): the branch reactance
            b (float): the branch susceptance
            status (bool): indicates if the branch is online or offline
            rateA (float): The 1st rating of the line.
            rateB (float): The 2nd rating of the line
            rateC (float): The 3rd rating of the line.
        """
        self.id = self._ids.__next__()

        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.
        
        self.from_bus = from_bus
        self.to_bus = to_bus
        self.r = r
        self.x = x
        self.b = 1/self.x
        self.g = 1/self.r
        self.G = self.r/(self.r**2 + self.x**2)
        self.B = -self.x/(self.r**2 + self.x**2)
        self.status = status
        self.rateA = rateA
        self.rateB = rateB
        self.rateC = rateC
        
        
    def diagonal_stamp(self, Y, J, from_node, to_node):
        Y[from_node][from_node] += self.G
        Y[from_node][to_node] += -self.B
        Y[to_node][from_node] += self.B
        Y[to_node][to_node] += self.G
        return Y, J
    
    def off_diagonal_stamp(self,Y,J,a,b,c,d):
        # nodes a,b,c,d make a square as so:
        # [a][c] | [a][d]
        # ---------------
        # [b][c] | [b][d]
        Y[a][c] += -self.G
        Y[a][d] += self.B
        Y[b][c] += -self.B
        Y[b][d] += -self.G
        return Y, J
    
    def stamp(self, Y, J):
        # diagonal stamp blocks
        Y, J = self.diagonal_stamp(Y,J,Buses.bus_map[self.from_bus].node_Vr,
                                   Buses.bus_map[self.from_bus].node_Vi)
        
        Y, J = self.diagonal_stamp(Y,J,Buses.bus_map[self.to_bus].node_Vr,
                                   Buses.bus_map[self.to_bus].node_Vi)
        # off-diagonal stamp blocks
        # top right
        Y, J = self.off_diagonal_stamp(Y,J,Buses.bus_map[self.from_bus].node_Vr,
                                   Buses.bus_map[self.from_bus].node_Vi,
                                   Buses.bus_map[self.to_bus].node_Vr,
                                   Buses.bus_map[self.to_bus].node_Vi)
        # bottom left
        Y, J = self.off_diagonal_stamp(Y,J,Buses.bus_map[self.to_bus].node_Vr,
                                   Buses.bus_map[self.to_bus].node_Vi,
                                   Buses.bus_map[self.from_bus].node_Vr,
                                   Buses.bus_map[self.from_bus].node_Vi)    
        
        return Y,J    
        
        