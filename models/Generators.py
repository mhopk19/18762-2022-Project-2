from __future__ import division
from itertools import count
from models.Buses import Buses
import scripts.global_vars as gv


class Generators:
    _ids = count(0)
    RemoteBusGens = dict()
    RemoteBusRMPCT = dict()
    gen_bus_key_ = {}
    total_P = 0
    base = gv.global_vars.MVAbase

    def __init__(self,
                 Bus,
                 P,
                 Vset,
                 Qmax,
                 Qmin,
                 Pmax,
                 Pmin,
                 Qinit,
                 RemoteBus,
                 RMPCT,
                 gen_type):
        """Initialize an instance of a generator in the power grid.

        Args:
            Bus (int): the bus number where the generator is located.
            P (float): the current amount of active power the generator is providing. [MW]
            Vset (float): the voltage setpoint that the generator must remain fixed at.
            Qmax (float): maximum reactive power [Mvar]
            Qmin (float): minimum reactive power [Mvar]
            Pmax (float): maximum active power [MW]
            
            Pmin (float): minimum active power [MW]
            Qinit (float): the initial amount of reactive power that the generator is supplying or absorbing.
            RemoteBus (int): the remote bus that the generator is controlling
            RMPCT (float): the percent of total MVAR required to hand the voltage at the controlled bus
            gen_type (str): the type of generator
        """

        self.id = self._ids.__next__()

        print("Bus:{}P:{}Vset:{}Qmax:{}Qmin:{}Pmax:{}\n Pmin:{}Qinit:{}RemoteBus:{}RMPCT:{}gen_type:{}".format(
            Bus, P, Vset, Qmax, Qmin, Pmax, Pmin, Qinit, RemoteBus, RMPCT, gen_type))
        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.
        self.Bus = Bus
        self.P = P / Generators.base
        self.Vset = Vset
        self.Qmax = Qmax / Generators.base
        self.Qmin = Qmin / Generators.base
        self.Pmax = Pmax / Generators.base
        self.Pmin = Pmin / Generators.base
        self.Qinit = Qinit / Generators.base
        self.RemoteBus = RemoteBus
        self.RMPCT = RMPCT
        self.gen_type = gen_type
        
    def dIrg_dVrg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        term1 = self.P/(Vrg**2 + Vig**2)
        term2 = -(2*Vrg*(self.P*Vrg - Q*Vig))/(Vrg**2 + Vig**2)**2
        return term1 + term2
    
    def dIrg_dVig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        term1 = Q/(Vrg**2 + Vig**2)
        term2 = -(2*Vig*(self.P*Vrg - Q*Vig))/(Vrg**2 + Vig**2)**2
        return term1 + term2
    
    def dIig_dVrg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        term1 = Q/(Vrg**2 + Vig**2)
        term2 = (2*Vrg*(self.P*Vig - Q*Vrg))/(Vrg**2 + Vig**2)**2
        return term1 + term2
    
    def dIig_dVig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        term1 = -self.P/(Vrg**2 + Vig**2)
        term2 = (2*Vig*(self.P*Vig - Q*Vrg))/(Vrg**2 + Vig**2)**2
        return term1 + term2
    
    def dIrg_dQg(self,Vrg,Vig):
        assert (Vrg!=0 or Vig!=0)
        num = -Vig
        denom = (Vig**2 + Vrg**2)
        return num/denom
    
    def dIig_dQg(self,Vrg,Vig):
        assert (Vrg!=0 or Vig!=0)
        num = Vrg
        denom = (Vig**2 + Vrg**2)
        return num/denom
    
    def Irg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        num = -self.P*Vrg - Q*Vig
        denom = (Vrg**2 + Vig**2)
        return num/denom
    
    def Iig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        num = -self.P*Vig + Q*Vrg
        denom = (Vrg**2 + Vig**2)
        return num/denom
    
    def stamp(self, Y, J, prev_v):
        # prev_v = Vrg_k,Vig_k,Qg_k
        #
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        q_node = Buses.bus_map[self.Bus].node_Q
        
        # conductance and VCVS
        # Rg and Ig differentials
        Y[v_node_r][v_node_r] += self.dIrg_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        Y[v_node_r][v_node_i] += self.dIrg_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        Y[v_node_i][v_node_i] += self.dIig_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        Y[v_node_i][v_node_r] += self.dIig_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        
        # Qg differentials
        Y[v_node_r][q_node] += self.dIrg_dQg(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i][q_node] += self.dIig_dQg(prev_v[v_node_r],prev_v[v_node_i])
        
        # historical values
        # Vrl
        J[v_node_r] -= self.Irg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) - \
            self.dIrg_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_r] -\
            self.dIrg_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_i] -\
            self.dIrg_dQg(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[q_node]
            
        # Vil
        J[v_node_r] -= self.Iig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) - \
            self.dIig_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_r] -\
            self.dIig_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_i] -\
            self.dIig_dQg(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[q_node]
        
        # V set row 
        # used to be multiples of 2 below
        v_eq_hist = -self.Vset**2 - prev_v[v_node_r]**2 - prev_v[v_node_i]**2   
        Y[q_node][v_node_r] += 2 * prev_v[v_node_r]
        Y[q_node][v_node_i] += 2 * prev_v[v_node_i]
        J[q_node] -= v_eq_hist 
        
            
        return Y, J
    