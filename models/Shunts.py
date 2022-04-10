from __future__ import division
from itertools import count
from models.Buses import Buses
import scripts.global_vars as gv

class Shunts:
    _ids = count(0)
    base = gv.global_vars.MVAbase
    
    def __init__(self,
                 Bus,
                 G_MW,
                 B_MVAR,
                 shunt_type,
                 Vhi,
                 Vlo,
                 Bmax,
                 Bmin,
                 Binit,
                 controlBus,
                 flag_control_shunt_bus=False,
                 Nsteps=[0],
                 Bstep=[0]):

        """ Initialize a shunt in the power grid.
        Args:
            Bus (int): the bus where the shunt is located
            G_MW (float): the active component of the shunt admittance as MW per unit voltage
            B_MVAR (float): reactive component of the shunt admittance as  MVar per unit voltage
            shunt_type (int): the shunt control mode, if switched shunt
            Vhi (float): if switched shunt, the upper voltage limit
            Vlo (float): if switched shunt, the lower voltage limit
            Bmax (float): the maximum shunt susceptance possible if it is a switched shunt
            Bmin (float): the minimum shunt susceptance possible if it is a switched shunt
            Binit (float): the initial switched shunt susceptance
            controlBus (int): the bus that the shunt controls if applicable
            flag_control_shunt_bus (bool): flag that indicates if the shunt should be controlling another bus
            Nsteps (list): the number of steps by which the switched shunt should adjust itself
            Bstep (list): the admittance increase for each step in Nstep as MVar at unity voltage
        """
        self.id = self._ids.__next__()

        self.Bus = Bus
        self.G_MW = G_MW / Shunts.base
        self.B_MVAR = B_MVAR / Shunts.base
        self.shunt_type = shunt_type
        self.Vhi = Vhi
        self.Vlo = Vlo
        self.Bmax = Bmax
        self.Bmin = Bmin
        self.Binit = Binit
        self.controlBus = controlBus
        self.flag_control_shunt_bus = flag_control_shunt_bus=False
        
        # You will need to implement the remainder of the __init__ function yourself.
        # You should also add some other class functions you deem necessary for stamping,
        # initializing, and processing results.
        
    def stamp(self,Y,J):
        v_node_r = Buses.bus_map[self.Bus].node_Vr
        v_node_i = Buses.bus_map[self.Bus].node_Vi
        
        # independent voltage source stamping
        Y[v_node_r][v_node_r] = self.G_MW
        Y[v_node_r][v_node_i] = -self.B_MVAR
        Y[v_node_i][v_node_r] = self.B_MVAR
        Y[v_node_i][v_node_i] = self.G_MW
        
        return Y,J  