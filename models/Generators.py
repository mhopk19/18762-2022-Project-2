from __future__ import division
from itertools import count
from models.Buses import Buses
import scripts.global_vars as gv

import sympy
from sympy import Symbol

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
        self.P_orig = P / Generators.base
        self.P = -self.P_orig #P / Generators.base
        self.Vset = Vset
        self.Qmax = Qmax / Generators.base
        self.Qmin = Qmin / Generators.base
        self.Pmax = Pmax / Generators.base
        self.Pmin = Pmin / Generators.base
        self.Qinit = Qinit / Generators.base
        self.RemoteBus = RemoteBus
        self.RMPCT = RMPCT
        self.gen_type = gen_type
        
    def assign_indexes(self):
        self.node_Vr = Buses.bus_map[self.Bus].node_Vr
        self.node_Vi = Buses.bus_map[self.Bus].node_Vi
        self.node_Q = Buses.bus_map[self.Bus].node_Q
        
    def assign_lamda_indexes(self):
        self.lambda_Vr = Buses.bus_map[self.Bus].lambda_Vr
        self.lambda_Vi = Buses.bus_map[self.Bus].lambda_Vi
        self.lambda_Q = Buses.bus_map[self.Bus].lambda_Q
        
    def dIrg_dVrg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        #term1 = -self.P/(Vrg**2 + Vig**2)
        #term2 = -(2*Vrg*(self.P*Vrg - Q*Vig))/(Vrg**2 + Vig**2)**2
        #return term1 + term2
        num = (self.P*(Vig**2 - Vrg**2) - 2 * Q * Vrg * Vig)
        denom = (Vrg ** 2 + Vig ** 2) ** 2
        return num / denom
    
    
    def dIrg_dVig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        #term1 = Q/(Vrg**2 + Vig**2)
        #term2 = -(2*Vig*(self.P*Vrg - Q*Vig))/(Vrg**2 + Vig**2)**2
        #return term1 + term2
        num = (Q * (Vrg**2 - Vig**2) - 2 * self.P * Vrg * Vig)
        denom = (Vrg ** 2 + Vig ** 2) ** 2
        return num / denom
    
    def dIig_dVrg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        #term1 = Q/(Vrg**2 + Vig**2)
        #term2 = (2*Vrg*(self.P*Vig - Q*Vrg))/(Vrg**2 + Vig**2)**2
        #return term1 + term2
        num = (Q * (Vrg**2 - Vig**2) - 2 * self.P * Vrg * Vig)
        denom = (Vrg ** 2 + Vig ** 2) ** 2
        return num / denom
    
    def dIig_dVig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        #term1 = -self.P/(Vrg**2 + Vig**2)
        #term2 = (2*Vig*(self.P*Vig - Q*Vrg))/(Vrg**2 + Vig**2)**2
        #return term1 + term2
        num = -(self.P*(Vig**2 - Vrg**2) - 2 * Q * Vrg * Vig)
        denom = (Vrg ** 2 + Vig ** 2) ** 2
        return num / denom
    
    def dIrg_dQg(self,Vrg,Vig):
        assert (Vrg!=0 or Vig!=0)
        num = Vig
        denom = (Vig**2 + Vrg**2)
        return num/denom
    
    def dIig_dQg(self,Vrg,Vig):
        assert (Vrg!=0 or Vig!=0)
        num = -Vrg
        denom = (Vig**2 + Vrg**2)
        return num/denom
    
    def Irg(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        num = self.P*Vrg + Q*Vig
        denom = (Vrg**2 + Vig**2)
        return num/denom
    
    def Iig(self,Vrg,Vig,Q):
        assert (Vrg!=0 or Vig!=0)
        num = self.P*Vig - Q*Vrg
        denom = (Vrg**2 + Vig**2)
        return num/denom
    
    def Q_constraint(self, Vrg, Vig):
        return self.Vset**2 - Vrg**2 - Vig**2
    
    def IIrg(self, Vrg, Vig, Q, L_r, L_i, L_q):
        LpQ = L_r * (self.Irg(Vrg,Vig,Q)) + L_i * (self.Iig(Vrg,Vig,Q)) + L_q * (self.Q_constraint(Vrg, Vig))
        return sympy.diff(LpQ, Vrg)
        
    def IIig(self, Vrg, Vig, Q, L_r, L_i, L_q):
        LpQ = L_r * (self.Irg(Vrg,Vig,Q)) + L_i * (self.Iig(Vrg,Vig,Q)) + L_q * (self.Q_constraint(Vrg, Vig))
        return sympy.diff(LpQ, Vig)
    
    def IIqg(self, Vrg, Vig, Q, L_r, L_i, L_q):
        LpQ = L_r * (self.Irg(Vrg,Vig,Q)) + L_i * (self.Iig(Vrg,Vig,Q)) + L_q * (self.Q_constraint(Vrg, Vig))
        return sympy.diff(LpQ, Q)
        
    
    def stamp(self, Y, J, prev_v):
        # prev_v = Vrg_k,Vig_k,Qg_k
        #
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        q_node = Buses.bus_map[self.Bus].node_Q
        
        # conductance and VCVS
        # Rg and Ig differentials
        Y[v_node_r ][v_node_r ] += self.dIrg_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        Y[v_node_r ][v_node_i ] += self.dIrg_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        
        Y[v_node_i ][v_node_i ] += self.dIig_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        Y[v_node_i ][v_node_r ] += self.dIig_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node])
        
        # Qg differentials
        Y[v_node_r ][q_node ] += self.dIrg_dQg(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i ][q_node ] += self.dIig_dQg(prev_v[v_node_r],prev_v[v_node_i])
        
        # historical values
        # Vrl
        J[v_node_r ] += -self.Irg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) + \
            self.dIrg_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_r] +\
            self.dIrg_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_i] +\
            self.dIrg_dQg(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[q_node]
            
        # Vil
        J[v_node_i ] += -self.Iig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) + \
            self.dIig_dVrg(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_r] +\
            self.dIig_dVig(prev_v[v_node_r],prev_v[v_node_i],prev_v[q_node]) * prev_v[v_node_i] +\
            self.dIig_dQg(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[q_node]
        
        # V set row 
        # used to be multiples of 2 below
        v_eq_hist = -self.Vset**2 - prev_v[v_node_r]**2 - prev_v[v_node_i]**2  
        #v_eq_hist = v_eq_hist -2 * prev_v[v_node_r]*prev_v[v_node_r] - 2 * prev_v[v_node_i]*prev_v[v_node_i]
        Y[q_node ][v_node_r ] += -2 * prev_v[v_node_r]
        Y[q_node ][v_node_i ] += -2 * prev_v[v_node_i]
        J[q_node ] += v_eq_hist 
        
            
        return Y, J
    
    def get_lagrange(self, symb_v, symb_lambda):
        # branch goes from node k to node m
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        Pg = Symbol('Pg{}'.format(self.id))
        Qg = Symbol('Qg{}'.format(self.id))
        
        # expression is in terms of currents
        denominator = (symb_v[v_node_r]**2 + symb_v[v_node_i]**2)
        r_expr = symb_lambda[v_node_r] * (-Pg * symb_v[v_node_r] - Qg * symb_v[v_node_i]) / denominator 
        i_expr = symb_lambda[v_node_i] *(-Pg * symb_v[v_node_i] - Qg * symb_v[v_node_r]) / denominator 

        
        return r_expr + i_expr
    
    def update_symbols(self, dict, prev_v):
        q_node = Buses.bus_map[self.Bus].node_Q
        dict['Pg{}'.format(self.id)] = self.P
        dict['Qg{}'.format(self.id)] = prev_v[q_node]
        
        return dict
    
    
    def stamp_dual(self, Y, J, prev_sol, size_Y, print_terms = False):
        dict = {}
        x_ind_end = Buses.lambda_v_start
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        q_node = Buses.bus_map[self.Bus].node_Q
        l_node_r = Buses.bus_map[self.Bus].lambda_Vr 
        l_node_i = Buses.bus_map[self.Bus].lambda_Vi
        l_node_q = Buses.bus_map[self.Bus].lambda_Q
        
        print("generator dual stamp nodes vr:{} vi:{} q_node:{} \nlr:{} li:{} lq:{}".format(v_node_r, v_node_i,
                                                                     q_node, l_node_r, 
                                                                     l_node_i, l_node_q))
        print("stamps for Generator:{}".format(self.id))
        
        dict["v_r"] = (sympy.Symbol("v_r"), prev_sol[v_node_r])
        dict["v_i"] = (sympy.Symbol("v_i"), prev_sol[v_node_i])
        dict["v_q"] = (sympy.Symbol("v_q"), prev_sol[q_node])
        dict["l_r"] = (sympy.Symbol("l_r"), prev_sol[l_node_r])
        dict["l_i"] = (sympy.Symbol("l_i"), prev_sol[l_node_i])
        dict["l_q"] = (sympy.Symbol("l_q"), prev_sol[l_node_q])
        
        # historical value calculations
        prev_IIrg = self.IIrg(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0], 
                              dict["l_r"][0], dict["l_i"][0], dict["l_q"][0])
        
        prev_IIig = self.IIig(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0], 
                              dict["l_r"][0], dict["l_i"][0], dict["l_q"][0])

        prev_IIqg = self.IIqg(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0], 
                              dict["l_r"][0], dict["l_i"][0], dict["l_q"][0])
        
        Beta_rk = self.IIrg(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0],
                            dict["l_r"][0], dict["l_i"][0], dict["l_q"][0]) - \
                dict["l_r"][0] * sympy.diff(prev_IIrg, dict["l_r"][0]) - \
                dict["l_i"][0] * sympy.diff(prev_IIrg, dict["l_i"][0]) - \
                dict["l_q"][0] * sympy.diff(prev_IIrg, dict["l_q"][0]) - \
                dict["v_r"][0] * sympy.diff(prev_IIrg, dict["v_r"][0]) - \
                dict["v_i"][0] * sympy.diff(prev_IIrg, dict["v_i"][0]) - \
                dict["v_q"][0] * sympy.diff(prev_IIrg, dict["v_q"][0])
                
        Beta_ik = self.IIig(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0],
                            dict["l_r"][0], dict["l_i"][0], dict["l_q"][0]) - \
                dict["l_r"][0] * sympy.diff(prev_IIig, dict["l_r"][0]) - \
                dict["l_i"][0] * sympy.diff(prev_IIig, dict["l_i"][0]) - \
                dict["l_q"][0] * sympy.diff(prev_IIig, dict["l_q"][0]) - \
                dict["v_r"][0] * sympy.diff(prev_IIig, dict["v_r"][0]) - \
                dict["v_i"][0] * sympy.diff(prev_IIig, dict["v_i"][0]) - \
                dict["v_q"][0] * sympy.diff(prev_IIig, dict["v_q"][0])
        
        veq_hist = self.IIqg(dict["v_r"][0], dict["v_i"][0], dict["v_q"][0],
                             dict["l_r"][0], dict["l_i"][0], dict["l_q"][0]) - \
                dict["l_r"][0] * sympy.diff(prev_IIqg, dict["l_r"][0]) - \
                dict["l_i"][0] * sympy.diff(prev_IIqg, dict["l_i"][0]) - \
                dict["l_q"][0] * sympy.diff(prev_IIqg, dict["l_q"][0]) - \
                dict["v_r"][0] * sympy.diff(prev_IIqg, dict["v_r"][0]) - \
                dict["v_i"][0] * sympy.diff(prev_IIqg, dict["v_i"][0]) - \
                dict["v_q"][0] * sympy.diff(prev_IIqg, dict["v_q"][0])
        
        
        # Real Stamps
        
        # real
        lambda_rr_next = sympy.diff(prev_IIrg, dict["l_r"][0])
        lambda_ri_next = sympy.diff(prev_IIrg, dict["l_i"][0])
        lambda_rq_next = sympy.diff(prev_IIrg, dict["l_q"][0])
        v_rr_next = sympy.diff(prev_IIrg, dict["v_r"][0])
        v_ri_next = sympy.diff(prev_IIrg, dict["v_i"][0])
        v_rq_next = sympy.diff(prev_IIrg, dict["v_q"][0])
        
        # imaginary
        lambda_ir_next = sympy.diff(prev_IIig, dict["l_r"][0])
        lambda_ii_next = sympy.diff(prev_IIig, dict["l_i"][0])
        lambda_iq_next = sympy.diff(prev_IIig, dict["l_i"][0])
        v_ir_next = sympy.diff(prev_IIig, dict["v_r"][0])
        v_ii_next = sympy.diff(prev_IIig, dict["v_i"][0])
        v_iq_next = sympy.diff(prev_IIig, dict["v_q"][0])
        
        # q-node
        lambda_qr_next = sympy.diff(prev_IIqg, dict["l_r"][0])
        lambda_qi_next = sympy.diff(prev_IIqg, dict["l_i"][0])
        lambda_qq_next = sympy.diff(prev_IIqg, dict["l_i"][0])
        v_qr_next = sympy.diff(prev_IIqg, dict["v_r"][0])
        v_qi_next = sympy.diff(prev_IIqg, dict["v_i"][0])
        v_qq_next = sympy.diff(prev_IIqg, dict["v_q"][0])  
        
              
        for i,item in dict.items():
            Beta_rk = Beta_rk.subs(item[0], item[1])
            Beta_ik = Beta_ik.subs(item[0], item[1])
            veq_hist = veq_hist.subs(item[0], item[1])
            
            lambda_rr_next = lambda_rr_next.subs(item[0], item[1])
            lambda_ri_next = lambda_ri_next.subs(item[0], item[1])
            lambda_rq_next = lambda_rq_next.subs(item[0], item[1])
            v_rr_next = v_rr_next.subs(item[0], item[1])
            v_ri_next = v_ri_next.subs(item[0], item[1])
            v_rq_next = v_rq_next.subs(item[0], item[1])

            lambda_ir_next = lambda_ir_next.subs(item[0], item[1])
            lambda_ii_next = lambda_ii_next.subs(item[0], item[1])
            lambda_iq_next = lambda_iq_next.subs(item[0], item[1])
            v_ir_next = v_ir_next.subs(item[0], item[1])
            v_ii_next = v_ii_next.subs(item[0], item[1])
            v_iq_next = v_iq_next.subs(item[0], item[1])

            lambda_qr_next = lambda_qr_next.subs(item[0], item[1])
            lambda_qi_next = lambda_qi_next.subs(item[0], item[1])
            lambda_qq_next = lambda_qq_next.subs(item[0], item[1])
            v_qr_next = v_qr_next.subs(item[0], item[1])
            v_qi_next = v_qi_next.subs(item[0], item[1])
            v_qq_next = v_qq_next.subs(item[0], item[1])            
            
            
        if (print_terms):
            print("Beta_rk", Beta_rk) 
            print("Beta_ik", Beta_ik)
            print("stamps for r row:l_r term {}\nl_i term {}\nv_r term {}\nv_i term {}".format(lambda_rr_next,
                                                                                    lambda_ri_next,
                                                   v_rr_next, v_ri_next))
            print("stamps for i row:l_r term {}\nl_i term {}\nv_r term {}\nv_i term {}".format(lambda_ir_next, 
                                                                                           lambda_ii_next,
                                                   v_ir_next, v_ii_next))
            
        Y[v_node_r][l_node_r] += lambda_rr_next
        Y[v_node_r][l_node_i] += lambda_ri_next
        Y[v_node_r][l_node_q] += lambda_rq_next
        Y[v_node_r][v_node_r] += v_rr_next
        Y[v_node_r][v_node_i] += v_ri_next
        Y[v_node_r][q_node] += v_rq_next
        
        Y[v_node_i][l_node_r] += lambda_ir_next
        Y[v_node_i][l_node_i] += lambda_ii_next
        Y[v_node_i][l_node_q] += lambda_iq_next
        Y[v_node_i][v_node_r] += v_ir_next
        Y[v_node_i][v_node_i] += v_ii_next
        Y[v_node_i][q_node] += v_rq_next
        
        Y[q_node][l_node_r] += lambda_qr_next
        Y[q_node][l_node_i] += lambda_qi_next
        Y[q_node][l_node_q] += lambda_qq_next
        Y[q_node][v_node_r] += v_qr_next
        Y[q_node][v_node_i] += v_qi_next
        Y[q_node][q_node] += v_qq_next
        
        J[v_node_r] += Beta_rk
        J[v_node_i] += Beta_ik
        J[q_node] += veq_hist

        return Y, J