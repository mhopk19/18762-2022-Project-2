import numpy as np
from models.Buses import Buses
import scipy
import scipy.sparse as sp
import scripts.sparse_matrices as spm
import scripts.sparse_matrices as sm

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
        self.generator_dict = {} # mapping generators to solution indices (generator, indices vr,vi,q)
        self.load_dict = {} # mapping load to solution indices
        self.slack_dict = {} # mapping slack to solution indices

        # need a log of the partials for calculations
        # load: real: Irl/Vrl, Irl/Vil PV bus imag: Iil/Vrl, Iil/Vil
        # PV bus:        
        self.partial_log = None
        self.solution_v = None
        
        # voltage limiting
        self.v_max_limit = 2
        self.v_min_limit = 0.0000000001
        self.delta_limit = 0.1
        
    def solve(self, Y, J, init_v):
        if (self.sparse == True):
            Y.generate_matrix_from_sparse()
            J.generate_matrix_from_sparse()
            init_v.generate_matrix_from_sparse()
            
            v_new = init_v.sparse_matrix - sp.linalg.inv(Y.sparse_matrix) @  \
                    (Y.sparse_matrix @ init_v.sparse_matrix - J.sparse_matrix)
            # fit the matrix output of todense() to an 1-D array format
            v_new = np.array(v_new.todense()).ravel()
            #print("v new:", v_new,"size",v_new.shape)
            # need to fit the matrix output into an array format
            return sm.sparse_vector(arr = v_new)
        else:
            rounded_y = np.round(np.matrix(Y).tolist(),2)
            rounded_y = Y
            #print('\n'.join([''.join(['{:8}'.format(item) for item in row]) 
            #                 for row in rounded_y]))
            
            v_new = init_v - np.linalg.inv(Y) @ (Y @ init_v - J)
        # calculate information for determining residuals
            return v_new

    def apply_limiting(self, v_sol, prev_v_sol, voltage_devices):
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
        for s in slack:
            self.Y, self.J = s.stamp(self.Y, self.J)
        for b in branch:
            self.Y, self.J = b.stamp(self.Y, self.J)
        for t in transformer:
            self.Y, self.J = t.stamp(self.Y, self.J)
        for sh in shunt:
            self.Y, self.J = sh.stamp(self.Y, self.J)
            

    def stamp_nonlinear(self, load, generator, prev_v):
        # loads
        for l in load:
            self.Y, self.J = l.stamp(self.Y, self.J, prev_v)
        # generators
        for g in generator:
            self.Y, self.J = g.stamp(self.Y, self.J, prev_v)
    
    def reset_stamps(self, size):
        if (self.sparse):
            self.Y = spm.sparse_matrix(size = size)
            self.J = spm.sparse_vector(size = size)
        else:
            self.Y = np.zeros((size,size))
            self.J = np.zeros(size)
        
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
        # TODO: PART 1, STEP 2.1 - Complete the stamp_linear function which stamps all linear power grid elements.
        #  This function should call the stamp_linear function of each linear element and return an updated Y matrix.
        #  You need to decide the input arguments and return values.
        self.stamp_linear(slack,branch,transformer,shunt)
        
        # linear stamps which can be used as the starting point in NR iterations
        linear_stamps = (self.Y, self.J)
        
        # # # Initialize While Loop (NR) Variables # # #
        # TODO: PART 1, STEP 2.2 - Initialize the NR variables
        err_max = np.inf  # maximum error at the current NR iteration
        tol = self.tol # chosen NR tolerance
        NR_count = 0  # current NR iteration
        
        # resizing sparse matrix to avoid singularity
        #if (sparse):
        #    self.Y.sparse_matrix.resize((self.size_y, self.size_y))
        #    self.J.sparse_matrix.resize((self.size_y, self.size_y))

        # # # Begin Solving Via NR # # #
        # TODO: PART 1, STEP 2.3 - Complete the NR While Loop
        while (err_max > tol and (NR_count < self.max_iters)):
            print("NR iteration: {}".format(NR_count))
            self.Y, self.J = linear_stamps
            
            # # # Stamp Nonlinear Power Grid Elements into Y matrix # # #
            # TODO: PART 1, STEP 2.4 - Complete the stamp_nonlinear function which stamps all nonlinear power grid
            #  elements. This function should call the stamp_nonlinear function of each nonlinear element and return
            #  an updated Y matrix. You need to decide the input arguments and return values.
            self.stamp_nonlinear(load, generator, v_sol)

            # # # Solve The System # # #
            # TODO: PART 1, STEP 2.5 - Complete the solve function which solves system of equations Yv = J. The
            #  function should return a new v_sol.
            #  You need to decide the input arguments and return values.
            prev_v_sol = v_sol
            v_sol = self.solve(self.Y, self.J, v_sol)

            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 1, STEP 2.6 - Finish the check_error function which calculates the maximum error, err_max
            #  You need to decide the input arguments and return values
            err_max = self.check_error(self.Y, self.J, v_sol)
            print("max error at iteration:{}".format(err_max))
            #print("solution vector: {}".format(v_sol))

            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 2, STEP 1 - Develop the apply_limiting function which implements voltage and reactive power
            #  limiting. Also, complete the else condition. Do not complete this step until you've finished Part 1.
            #  You need to decide the input arguments and return values.
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

        return v_sol
