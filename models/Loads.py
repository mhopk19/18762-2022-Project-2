from __future__ import division
from itertools import count
from models.Buses import Buses
import math
import scripts.global_vars as gv

class Loads:
    _ids = count(0)
    base = gv.global_vars.MVAbase

    def __init__(self,
                 Bus,
                 P,
                 Q,
                 IP,
                 IQ,
                 ZP,
                 ZQ,
                 area,
                 status):
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
        self.id = Loads._ids.__next__()

        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.

        self.Bus = Bus
        self.P = P / Loads.base
        self.Q = Q / Loads.base
        self.IP = IP
        self.IQ = IQ
        self.ZP = ZP
        self.ZQ = ZQ
        self.area = area
        self.status = status
        
    def dIrl_dVrl(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        term1 = self.P/(Vrl**2 + Vil**2)
        term2 = -(2*Vrl*(self.P*Vrl + self.Q*Vil))/(Vrl**2 + Vil**2)**2
        return term1 + term2
    
    def dIrl_dVil(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        term1 = self.Q/(Vrl**2 + Vil**2)
        term2 = -(2*Vil*(self.P*Vrl + self.Q*Vil))/(Vrl**2 + Vil**2)**2
        return term1 + term2
    
    def dIil_dVrl(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        term1 = -self.Q/(Vrl**2 + Vil**2)
        term2 = -(2*Vrl*(self.P*Vil - self.Q*Vrl))/(Vrl**2 + Vil**2)**2
        return term1 + term2
    
    def dIil_dVil(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        term1 = self.P/(Vrl**2 + Vil**2)
        term2 = -(2*Vil*(self.P*Vil - self.Q*Vrl))/(Vrl**2 + Vil**2)**2
        return term1 + term2
    
    def Irl(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        num = self.P*Vrl + self.Q*Vil
        denom = (Vrl**2 + Vil**2)
        return num/denom
    
    def Iil(self,Vrl,Vil):
        assert (Vrl!=0 or Vil!=0)
        num = self.P*Vil - self.Q*Vrl
        denom = (Vrl**2 + Vil**2)
        return num/denom
    
    def stamp(self, Y, J, prev_v):
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        # conductance and VCVS
        Y[v_node_r][v_node_r] += self.dIrl_dVrl(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_r][v_node_i] += self.dIrl_dVil(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i][v_node_r] += self.dIil_dVrl(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i][v_node_i] += self.dIil_dVil(prev_v[v_node_r],prev_v[v_node_i])
        
        # historical values
        # Vrl
        J[v_node_r] -= self.Irl(prev_v[v_node_r],prev_v[v_node_i]) - \
            self.dIrl_dVrl(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_r] -\
            self.dIrl_dVil(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_i]
            
        # Vil
        J[v_node_i] -= self.Iil(prev_v[v_node_r],prev_v[v_node_i]) - \
            self.dIil_dVrl(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_r] -\
            self.dIil_dVil(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_i]
            
        return Y, J
        
        
        
        
    