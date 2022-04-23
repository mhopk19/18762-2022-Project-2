from __future__ import division
from itertools import count
from models.Buses import Buses


class Transformers:
    _ids = count(0)

    def __init__(self,
                 from_bus,
                 to_bus,
                 r,
                 x,
                 status,
                 tr,
                 ang,
                 Gsh_raw,
                 Bsh_raw,
                 rating):
        """Initialize a transformer instance

        Args:
            from_bus (int): the primary or sending end bus of the transformer.
            to_bus (int): the secondary or receiving end bus of the transformer
            r (float): the line resitance of the transformer in
            x (float): the line reactance of the transformer
            status (int): indicates if the transformer is active or not
            tr (float): transformer turns ratio
            ang (float): the phase shift angle of the transformer
            Gsh_raw (float): the shunt conductance of the transformer
            Bsh_raw (float): the shunt admittance of the transformer
            rating (float): the rating in MVA of the transformer
        """
        self.id = self._ids.__next__()

        self.from_bus = from_bus
        self.to_bus = to_bus
        self.r = r
        self.x = x
        self.status = status
        self.tr = tr
        self.ang = ang
        self.Gsh_raw = Gsh_raw
        self.Bsh_raw = Bsh_raw
        self.rating = rating
        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.
    
    def diagonal_stamp(self, Y, J, from_node, to_node):
        Y[from_node][from_node] += self.tr * self.x + (self.tr * (self.tr - 1)) * self.x
        Y[from_node][to_node] += -self.tr * self.x
        Y[to_node][from_node] += -self.tr * self.x
        Y[to_node][to_node] += self.tr * self.x + (1 - self.tr) * self.x
        return Y, J
    
    def stamp(self, Y, J):
        """
        Y, J = self.diagonal_stamp(Y,J,Buses.bus_map[self.from_bus].node_Vr,
                                   Buses.bus_map[self.to_bus].node_Vr) 
        
        Y, J = self.diagonal_stamp(Y,J,Buses.bus_map[self.from_bus].node_Vi,
                                   Buses.bus_map[self.to_bus].node_Vi)
        """
        
        # diagonal stamp blocks
        kr = Buses.bus_map[self.from_bus].node_Vr
        ki = Buses.bus_map[self.from_bus].node_Vi
        mr = Buses.bus_map[self.to_bus].node_Vr
        mi = Buses.bus_map[self.to_bus].node_Vi
        
        c = 1/(self.r**2+self.x**2) # constant which appears in each term
        
        # terms are coming from I_km,r equation
        Y[kr ][kr ] += c*self.r*self.tr**2
        Y[kr ][ki ] += c*self.x*self.tr**2
        Y[kr ][mr ] += -c*self.tr*self.r
        Y[kr ][mi ] += -c*self.tr*self.x
        # terms are coming I_km,i equation
        Y[ki ][ki ] += c*self.r*self.tr**2
        Y[ki ][kr ] += -c*self.x*self.tr**2
        Y[ki ][mr ] += c*self.tr*self.x
        Y[kr ][mi ] += -c*self.tr*self.r   
        
        # terms are coming from I_mk,r equation
        Y[mr ][kr ] += c*self.tr*self.r
        Y[mr ][ki ] += c*self.tr*self.x
        Y[mr ][mr ] += c*self.r
        Y[mr ][mi ] += c*self.x
        
        # terms are coming from I_mk,i equation
        Y[mi ][kr ] += -c*self.tr*self.x
        Y[mi ][ki ] += c*self.tr*self.r
        Y[mi ][mr ] += -c*self.x
        Y[mi ][mi ] += c*self.r  
        
        return Y,J  
        
        