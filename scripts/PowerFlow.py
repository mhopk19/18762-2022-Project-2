import numpy as np
from models.Buses import Buses
import scipy
import scipy.sparse as sp
import scripts.sparse_matrices as spm
import scripts.sparse_matrices as sm

import sympy
from sympy import symbols

import matplotlib.pyplot as plt

def print_matrix(matrix, width = 4):
    matrix_list = np.round(matrix, width).tolist()
    for row in matrix_list:
        str = ""
        for item in row:
            if ((item > 0) or (item < 0)):
                # maintain a maxium amount of decimal places
                # account for negative signs
                neg_flag = (item < 0)
                if (neg_flag):
                    int_width = int(np.floor(np.log10(item*-1)))   
                else:
                    int_width = int(np.floor(np.log10(item)))
                    
                decimal_places = '{}'.format(max(0, width - neg_flag - int_width))
            else:
                decimal_places = '{}'.format(width)
                
            substr = '{:.'+decimal_places+'f} '
            str = str+substr.format(item)
        print(str)
    #print('\n'.join([''.join(['{:.5f}'.format(item) for item in row]) 
    #            for row in matrix_list]))


class PowerFlow:

    def __init__(self,
                 case_name,
                 tol,
                 max_iters,
                 enable_limiting,
                 sparse,
                 size_y):
        """Initialize the PowerFlow instance.

        Args:
            case_name (str): A string with the path to the test case.
            tol (float): The chosen NR tolerance.
            max_iters (int): The maximum number of NR iterations.
            enable_limiting (bool): A flag that indicates if we use voltage limiting or not in our solver.
        """
        # Clean up the case name string
        case_name = case_name.replace('.RAW', '')
        case_name = case_name.replace('testcases/', '')

        self.case_name = case_name
        self.tol = tol
        self.max_iters = max_iters
        self.enable_limiting = enable_limiting
        self.sparse = sparse
        
        # added variables
        self.Y = None
        self.J = None
        self.size_y = size_y
        
        self.bus_map = {} # tuple mapping from integer number to bus

        # need a log of the partials for calculations
        # load: real: Irl/Vrl, Irl/Vil PV bus imag: Iil/Vrl, Iil/Vil
        # PV bus:        
        self.partial_log = None
        self.solution_v = None
        
        # voltage limiting
        self.v_max_limit = 2
        self.v_min_limit = 0.0000000001
        self.delta_limit = 0.1
        
        # feasibility variables
        self.feasibility_dict = {}
        self.lgr = None
        self.dx_lgr = None
        self.symb_v = None
        self.symb_lambda = None
        self.symb_dict = {}
        self.H_grad = None # the sub matrix of Y used in feasibility analysis
        self.L_grad = None

           
    def solve(self, Y, J, init_v):
        if (self.sparse == True):
            Y.generate_matrix_from_sparse()
            J.generate_matrix_from_sparse()
            init_v.generate_matrix_from_sparse()
            # NR step
            v_new = init_v.sparse_matrix - sp.linalg.inv(Y.sparse_matrix) @  \
                    (Y.sparse_matrix @ init_v.sparse_matrix - J.sparse_matrix)
            # fit the matrix output of todense() to an 1-D array format
            v_new = np.array(v_new.todense()).ravel()
            
            return sm.sparse_vector(arr = v_new)
        else:
            rounded_y = np.round(np.matrix(Y).tolist(),2)
            rounded_y = Y
            
            v_new = init_v - np.linalg.inv(Y) @ (Y @ init_v - J)
        # calculate information for determining residuals
            return v_new

    def apply_limiting(self, v_sol, prev_v_sol, voltage_devices):
        if (self.sparse):
            v_sol = np.array(v_sol.todense()).ravel()
            prev_v_sol = np.array(prev_v_sol.todense()).ravel()
            print("vsol:{} pre_v_sol:{}".format(v_sol, prev_v_sol))
        
        limit_vector = np.array(np.size(v_sol) * [self.delta_limit])
        # calculate delta v as well as its magnitude and sign
        delta_v = v_sol - prev_v_sol
        sign_delta_v = np.sign(delta_v)
        norm_delta_v = np.abs(delta_v)
        # calculate the voltage limited value of  v
        new_v = (sign_delta_v * np.minimum(norm_delta_v,limit_vector)) + prev_v_sol
        
        # perform absolute limits on all values
        new_v = np.minimum(new_v, np.size(v_sol) * [self.v_max_limit] )
        new_v = np.maximum(new_v, np.size(v_sol) * [self.v_min_limit] )
        
        # for all values that are NOT voltages replace new_delta_v values with v_sol values
        voltage_indices = {}
        # collect all the indices of v_sol that are voltages
        for device in voltage_devices:
            voltage_indices[Buses.bus_map[device.Bus].node_Vr] = 1
            voltage_indices[Buses.bus_map[device.Bus].node_Vi] = 1
        # subtract sets to find the indices that do not correspond to voltages
        voltage_indices = set(voltage_indices)
        all_indices = set(np.arange(0,np.size(v_sol),1))
        non_voltage_indices = all_indices.difference(voltage_indices)
        #print("voltage indices", voltage_indices)
        #print("all indices", all_indices)
        #print("non voltage indices", non_voltage_indices)
        
        # recover old values for the non_voltage indices
        for ind in non_voltage_indices:
            new_v[ind] = v_sol[ind]
        
        print("new v", new_v)
        
        if (self.sparse):
            return sm.sparse_vector(arr =new_v)
        else:
            return new_v

    def check_error(self, Y, J, v_sol):
        if (self.sparse == True):
            Y.generate_matrix_from_sparse()
            J.generate_matrix_from_sparse()
            v_sol.generate_matrix_from_sparse()
                      
            err_vector = (Y.sparse_matrix @ v_sol.sparse_matrix) - J.sparse_matrix
            err_max = np.max(np.abs(np.array(err_vector.todense()).ravel()))
            return err_max
        else:
            #print("Y",Y)
            #print("v sol", v_sol)
            #print("YV", Y @ v_sol)
            #print("J", J)
            err_vector = (Y @ v_sol) - J
            err_max = np.max(err_vector)
            return err_max        

    def get_hist_vars(self):
        """
        returns historical variables V r,g hist, V i,g hist and V e,q hist
        using self.last_v a record of the last solution
        """
        pass
        
    def stamp_linear(self, slack, branch, transformer, shunt):
        for comp in slack + branch + transformer + shunt:
            self.Y, self.J = comp.stamp(self.Y, self.J)
            
    def stamp_nonlinear(self, load, generator, prev_v):
        for comp in load + generator:
            self.Y, self.J = comp.stamp(self.Y, self.J, prev_v)
    
    def reset_stamps(self, size):
        if (self.sparse):
            self.Y = spm.sparse_matrix(size = size)
            self.J = spm.sparse_vector(size = size)
        else:
            self.Y = np.zeros((size,size))
            self.J = np.zeros(size)
            
    def reset_stamps_feasibility(self,size_y,size_x):
        if (self.sparse):
            self.Y = spm.sparse_matrix(size = size_y)
            self.J = spm.sparse_vector(size = size_y)
        else:
            self.Y = np.zeros((size_y,size_y))
            self.H_grad = np.zeros((size_x, size_x))
            self.lower_J = np.zeros(size_x)
            self.L_grad = np.zeros((size_y - size_x, size_y - size_x))
            self.J = np.zeros(size_y)  
                  
        
    def run_powerflow(self,
                      v_init,
                      bus,
                      slack,
                      generator,
                      transformer,
                      branch,
                      shunt,
                      load):
        """Runs a positive sequence power flow using the Equivalent Circuit Formulation.

        Args:
            v_init (np.array): The initial solution vector which has the same number of rows as the Y matrix.
            bus (list): Contains all the buses in the network as instances of the Buses class.
            slack (list): Contains all the slack generators in the network as instances of the Slack class.
            generator (list): Contains all the generators in the network as instances of the Generators class.
            transformer (list): Contains all the transformers in the network as instance of the Transformers class.
            branch (list): Contains all the branches in the network as instances of the Branches class.
            shunt (list): Contains all the shunts in the network as instances of the Shunts class.
            load (list): Contains all the loads in the network as instances of the Load class.

        Returns:
            v(np.array): The final solution vector.

        """
        # create bus mapping
        for b in bus:
            self.bus_map[b.Bus] = b

        step = 0
        
        # # # Copy v_init into the Solution Vectors used during NR, v, and the final solution vector v_sol # # #
        v = np.copy(v_init)
        v_sol = np.copy(v)
        if (self.sparse == True):
            v_sol = sm.sparse_vector(arr = v_init)
        v_size = self.size_y
        self.solution_v = np.copy(v)

        # initializing the MNA matrix/vectors
        self.reset_stamps(v_size)

        # # # Stamp Linear Power Grid Elements into Y matrix # # #
        for comp in slack + branch + transformer + shunt:
            self.Y, self.J = comp.stamp(self.Y, self.J)
        
        # linear stamps which can be used as the starting point in NR iterations
        linear_stamps = (self.Y, self.J)
        
        # # # Initialize While Loop (NR) Variables # # #
        err_max = np.inf  # maximum error at the current NR iteration
        tol = self.tol # chosen NR tolerance
        NR_count = 0  # current NR iteration
        
        # TODO: PART 1, STEP 2.3 - Complete the NR While Loop
        while (err_max > tol and (NR_count < self.max_iters)):
            print("NR iteration: {}".format(NR_count))
            self.Y, self.J = linear_stamps
            
            # # # Stamp Nonlinear Power Grid Elements into Y matrix # # #
            for comp in load + generator:
                self.Y, self.J = comp.stamp(self.Y, self.J, v_sol)

            # # # Solve The System # # #
            prev_v_sol = v_sol
            v_sol = self.solve(self.Y, self.J, v_sol)

            # # # Compute The Error at the current NR iteration # # #
            err_max = self.check_error(self.Y, self.J, v_sol)
            print("max error at iteration:{}".format(err_max))
            #print("solution vector: {}".format(v_sol))

            # # # Compute The Error at the current NR iteration # # #
            if self.enable_limiting and err_max > tol:
                print("enable limiting")
                v_sol = self.apply_limiting(v_sol, prev_v_sol, generator + slack + load)
            else:
                if (self.enable_limiting):
                    v_sol = self.apply_limiting(v_sol, prev_v_sol, generator + slack + load)
                    err_max = self.check_error(self.Y, self.J, v_sol)
            
            print("NR iteration", NR_count)
            prev_v_sol = v_sol            
            NR_count = NR_count + 1

        return v_sol, NR_count

    def run_feasibility_powerflow(self,
                      v_init,
                      bus,
                      slack,
                      generator,
                      transformer,
                      branch,
                      shunt,
                      load,
                      injection):
        
        def plot_matrix(matrix):
            plt.imshow(matrix)
            plt.show()
        
        n_bus = 0
        # create bus mapping
        for b in bus:
            n_bus = n_bus + 1
            self.bus_map[b.Bus] = b

        step = 0
        
        # # # Copy v_init into the Solution Vectors used during NR, v, and the final solution vector v_sol # # #
        v = np.copy(v_init)
        v_sol = np.copy(v)
        if (self.sparse == True):
            v_sol = sm.sparse_vector(arr = v_init)
        # actual v_size is different than size of the full matrix due to feasibility currents
        v_size = self.size_y 
        x_size = self.size_y - (2 * n_bus)
        lambda_size = (2 * n_bus)
        
        self.solution_v = np.copy(v)

        # initializing the MNA matrix/vectors
        self.reset_stamps_feasibility(v_size, x_size)
        
        # STAMPING LINEAR COMPONENTS
        
        # stamping cross-term Jacobian
        for comp in slack + branch + transformer + shunt:
            self.H_grad, self.lower_J = comp.stamp(self.H_grad, self.lower_J)
        
        # linear stamps which can be used as the starting point in NR iterations
        H_grad_linear = self.H_grad
        lower_J_linear = self.lower_J
        
        print("v_size: {} x_size:{} lamda_size:{}".format(v_size, x_size, lambda_size))
        print("Y sub1", self.Y[v_size - x_size:v_size , 0:x_size].shape)
        print("Y sub2", self.Y[0:x_size , v_size - x_size:v_size].shape)
        print("Hgrad", H_grad_linear.shape)
        
        # CREATE Y_linear
        
        # add cross term Jacobians
        # grad H
        self.Y[v_size - x_size:v_size , 0:x_size] = H_grad_linear
        self.J[v_size - x_size:v_size] = lower_J_linear
        
        
        # Lagrangian Hessian !!
        
        Y_linear, J_linear = (self.Y, self.J)
        
        print("linear Y", Y_linear)
        
        # CREATE J_linear
        
        
        # initialize feasibility analysis matrices
        
        
        plot_matrix(Y_linear)
    
        # # # Set Hyper-parameters
        tol = self.tol # chosen NR tolerance
        NR_count = 0  # current NR iteration
        err_max = np.inf # initial error to start loop
        while (err_max > tol and (NR_count < self.max_iters)):
            # # # ADDING THE NON-LINEAR COMPONENTS
            
            # # stamping H component of Y matrix
            for comp in load + generator:
                self.H_grad, self.lower_J = comp.stamp(H_grad_linear, self.lower_J, v_sol)
            for comp in injection:
                _, self.lower_J = comp.stamp(self.Y, self.lower_J, v_sol[v_size-x_size:v_size])
            
            # add H components to the Y matrix
            self.Y[v_size - x_size:v_size , 0:x_size] = self.H_grad
            self.J[v_size - x_size:v_size] = self.lower_J
            
            # # stamping upper part of Y matrix
            
            plot_matrix(self.Y)
            print_matrix(self.Y)
            print("J array", self.J)
            
            # # # Solve The System # # #
            prev_v_sol = v_sol
            #v_sol = self.solve(self.Y, self.J, v_sol)

            # # # Compute The Error at the current NR iteration # # #
            
            err_max = 0
            #err_max = self.check_error(self.Y, self.J, v_sol)
            
            #print("max error at iteration:{}".format(err_max))
            #print("solution vector: {}".format(v_sol))

            
            print("NR iteration: {}".format(NR_count))
            prev_v_sol = v_sol            
            NR_count = NR_count + 1
        
        return v_sol, NR_count


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