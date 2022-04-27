from __future__ import division
from itertools import count
from models.Buses import Buses
import math
import scripts.global_vars as gv

import sympy
from sympy import Symbol


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
    
    # dual expressions
    
    def IIrl(self,Vrl,Vil,L_r,L_i):
        LpQ = L_r * (self.Irl(Vrl,Vil)) + L_i * (self.Iil(Vrl,Vil))
        return sympy.diff(LpQ, Vrl)
        
    def IIil(self,Vrl,Vil,L_r,L_i):
        LpQ = L_r * (self.Irl(Vrl,Vil)) + L_i * (self.Iil(Vrl,Vil))
        return sympy.diff(LpQ, Vil)
        
    
    def stamp(self, Y, J, prev_v):
        
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        # conductance and VCVS
        Y[v_node_r ][v_node_r ] += self.dIrl_dVrl(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_r ][v_node_i ] += self.dIrl_dVil(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i ][v_node_r ] += self.dIil_dVrl(prev_v[v_node_r],prev_v[v_node_i])
        Y[v_node_i ][v_node_i ] += self.dIil_dVil(prev_v[v_node_r],prev_v[v_node_i])
        
        # historical values
        # Vrl
        J[v_node_r ] -= self.Irl(prev_v[v_node_r],prev_v[v_node_i]) - \
            self.dIrl_dVrl(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_r] -\
            self.dIrl_dVil(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_i]
            
        # Vil
        J[v_node_i ] -= self.Iil(prev_v[v_node_r],prev_v[v_node_i]) - \
            self.dIil_dVrl(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_r] -\
            self.dIil_dVil(prev_v[v_node_r],prev_v[v_node_i]) * prev_v[v_node_i]
            
        return Y, J
        
    def get_lagrange(self, symb_v, symb_lambda):
        # branch goes from node k to node m
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi

        Pl = Symbol('Pl{}'.format(self.id))
        Ql = Symbol('Ql{}'.format(self.id))
        
        # expression is in terms of currents
        denominator = (symb_v[v_node_r]**2 + symb_v[v_node_i]**2)
        r_expr = symb_lambda[v_node_r] * (Pl * symb_v[v_node_r] + Ql * symb_v[v_node_i]) / denominator 
        i_expr = symb_lambda[v_node_i] *(Pl * symb_v[v_node_i] - Ql * symb_v[v_node_r]) / denominator 


        
        return r_expr + i_expr       
        
    def update_symbols(self, dict, prev_v):
        dict['Pl{}'.format(self.id)] = self.P
        dict['Ql{}'.format(self.id)] = self.Q
        
        return dict
    
    def stamp_dual(self, Y, J, prev_sol, size_Y):
        dict = {}
        x_ind_end = size_Y - 2 * (self.Bus)
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        l_node_r = 2 * (self.Bus - 1) 
        l_node_i = 2 * (self.Bus - 1) + 1
        
        print("stamps for Load:{}".format(self.id))
        
        dict["v_r"] = (sympy.Symbol("v_r"), prev_sol[v_node_r])
        dict["v_i"] = (sympy.Symbol("v_i"), prev_sol[v_node_i])
        dict["l_r"] = (sympy.Symbol("l_r"), prev_sol[l_node_r])
        dict["l_i"] = (sympy.Symbol("l_i"), prev_sol[l_node_i])
        
        # historical value calculations
        prev_IIrl = self.IIrl(dict["v_r"][0], dict["v_i"][0], dict["l_r"][0], dict["l_i"][0])
        prev_IIil = self.IIil(dict["v_r"][0], dict["v_i"][0], dict["l_r"][0], dict["l_i"][0])
        
        Beta_rk = self.IIrl(dict["v_r"][0], dict["v_i"][0], dict["l_r"][0], dict["l_i"][0]) - \
                dict["l_r"][0] * sympy.diff(prev_IIrl, dict["l_r"][0]) - \
                dict["l_i"][0] * sympy.diff(prev_IIrl, dict["l_i"][0]) - \
                dict["v_r"][0] * sympy.diff(prev_IIrl, dict["v_r"][0]) - \
                dict["v_i"][0] * sympy.diff(prev_IIrl, dict["v_i"][0])
                
        Beta_ik = self.IIil(dict["v_r"][0], dict["v_i"][0], dict["l_r"][0], dict["l_i"][0]) - \
                dict["l_r"][0] * sympy.diff(prev_IIil, dict["l_r"][0]) - \
                dict["l_i"][0] * sympy.diff(prev_IIil, dict["l_i"][0]) - \
                dict["v_r"][0] * sympy.diff(prev_IIil, dict["v_r"][0]) - \
                dict["v_i"][0] * sympy.diff(prev_IIil, dict["v_i"][0])
        
        # Real Stamps
        
        # real
        lambda_rr_next = sympy.diff(prev_IIrl, dict["l_r"][0])
        lambda_ri_next = sympy.diff(prev_IIrl, dict["l_i"][0])
        v_rr_next = sympy.diff(prev_IIrl, dict["v_r"][0])
        v_ri_next = sympy.diff(prev_IIrl, dict["v_i"][0])
        
        # imaginary
        lambda_ir_next = sympy.diff(prev_IIil, dict["l_r"][0])
        lambda_ii_next = sympy.diff(prev_IIil, dict["l_i"][0])
        v_ir_next = sympy.diff(prev_IIil, dict["v_r"][0])
        v_ii_next = sympy.diff(prev_IIil, dict["v_i"][0])
        
        for i,item in dict.items():
            Beta_rk = Beta_rk.subs(item[0], item[1])
            Beta_ik = Beta_ik.subs(item[0], item[1])
            
            lambda_rr_next = lambda_rr_next.subs(item[0], item[1])
            lambda_ri_next = lambda_ri_next.subs(item[0], item[1])
            v_rr_next = v_rr_next.subs(item[0], item[1])
            v_ri_next = v_ri_next.subs(item[0], item[1])

            lambda_ir_next = lambda_ir_next.subs(item[0], item[1])
            lambda_ii_next = lambda_ii_next.subs(item[0], item[1])
            v_ir_next = v_ir_next.subs(item[0], item[1])
            v_ii_next = v_ii_next.subs(item[0], item[1])
            
        print("Beta_rk", Beta_rk) 
        print("Beta_ik", Beta_ik)
        print("stamps for r row:l_r term {}\nl_i term {}\nv_r term {}\nv_i term {}".format(lambda_rr_next,
                                                                                lambda_ri_next,
                                               v_rr_next, v_ri_next))
        print("stamps for i row:l_r term {}\nl_i term {}\nv_r term {}\nv_i term {}".format(lambda_ir_next, 
                                                                                       lambda_ii_next,
                                               v_ir_next, v_ii_next))
        
        Y[l_node_r][x_ind_end + l_node_r] += lambda_rr_next
        Y[l_node_r][x_ind_end + l_node_i] += lambda_ri_next
        Y[l_node_r][v_node_r] += v_rr_next
        Y[l_node_r][v_node_i] += v_ri_next
        
        Y[l_node_r][x_ind_end + l_node_r] += lambda_ir_next
        Y[l_node_r][x_ind_end + l_node_i] += lambda_ii_next
        Y[l_node_r][v_node_r] += v_ir_next
        Y[l_node_r][v_node_i] += v_ii_next
        
        J[l_node_r] += Beta_rk
        J[l_node_i] += Beta_ik

        return Y, J
         