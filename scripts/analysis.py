import sympy

def analyze_system(self, v_init, bus, slack, generator, transformer,
                   branch, shunt, load, injection):
    v = v_init
    v_size = len(v_init)
    lambda_size = 2 * len(bus)
    # feasibility calculations
    components = branch + generator + load
    self.symb_v = symbols('v:{}'.format(v_size))
    self.symb_lambda = symbols('l:{}'.format(lambda_size))
    self.Lgr = 0
    
    for comp in [generator[0]]:
        self.Lgr = self.Lgr + comp.get_lagrange(self.symb_v, self.symb_lambda)
    
    # update dictionary
    for comp in [generator[0]]:
        self.symb_dict = comp.update_symbols(self.symb_dict, v)
        
    print("Lagrangian: {}".format(self.Lgr))
    #print("After substitution: {}".format(self.Lgr.subs(self.symb_dict)))
    
    #print("Buses map", Buses.bus_map)
    #for _,i in Buses.bus_map.items():
    #    print("r",i.node_Vr)
    #    print("i",i.node_Vi)
    
    test_comp = generator[0]
    v_node_r = Buses.bus_map[test_comp.Bus].node_Vr
    v_node_i = Buses.bus_map[test_comp.Bus].node_Vi
    imag_deriv = sympy.diff(self.Lgr, "v{}".format(v_node_i))
    real_deriv = sympy.diff(self.Lgr, "v{}".format(v_node_r))
    
    print("load: vr index:{} vi index:{}".format(v_node_r, v_node_i))
    print("real derivative", imag_deriv)
    print("imaginary derivative", real_deriv)
    
    print("second imaginary derivative", sympy.diff(imag_deriv,"v{}".format(v_node_i)))
    print("second real derivative", sympy.diff(real_deriv,"v{}".format(v_node_r)))