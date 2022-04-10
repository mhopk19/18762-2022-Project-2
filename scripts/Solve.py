from parsers.parser import parse_raw
from scripts.PowerFlow import PowerFlow
from scripts.process_results import process_results
from scripts.initialize import initialize
from models.Buses import Buses
import scripts.sparse_matrices as sm
import numpy as np
import time

def solve(TESTCASE, SETTINGS):
    """Run the power flow solver.

    Args:
        TESTCASE (str): A string with the path to the test case.
        SETTINGS (dict): Contains all the solver settings in a dictionary.

    Returns:
        None
    """
    # TODO: PART 1, STEP 0 - Initialize all the model classes in the models directory (models/) and familiarize
    #  yourself with the parameters of each model. Use the docs/DataFormats.pdf for assistance.

    # # # Parse the Test Case Data # # #
    case_name = TESTCASE
    parsed_data = parse_raw(case_name)

    # # # Assign Parsed Data to Variables # # #
    bus = parsed_data['buses']
    slack = parsed_data['slack']
    generator = parsed_data['generators']
    transformer = parsed_data['xfmrs']
    branch = parsed_data['branches']
    shunt = parsed_data['shunts']
    load = parsed_data['loads']

    # # # Solver Settings # # #
    tol = SETTINGS['Tolerance']  # NR solver tolerance
    max_iters = SETTINGS['Max Iters']  # maximum NR iterations
    enable_limiting = SETTINGS['Limiting']  # enable/disable voltage and reactive power limiting
    sparse = SETTINGS['Sparse']

    # # # Assign System Nodes Bus by Bus # # #
    # We can use these nodes to have predetermined node number for every node in our Y matrix and J vector.
    for ele in bus:
        ele.assign_nodes()

    # Assign any slack nodes
    for ele in slack:
        ele.assign_nodes()
    
    # debugging
    print("bus:{}\nslack:{}\ngenerator:{}".format(bus,slack,generator))
    time.sleep(10)
    print("transformer:{}\nbranch:{}\nshunt:{}\nload:{}".format(transformer,branch,shunt,load))
    time.sleep(10)
    
    for g in generator:
        print("g:{}\n{}".format(g,vars(g)))

    for s in slack:
        print("s:{}\n{}".format(s,vars(s)))
        
    for l in load:
        print("l:{}\n{}".format(l,vars(l)))
        
    for b in bus:
        print("b:{}\n{}".format(b,vars(b)))     
        
    for br in branch:
        print("br:{}\n{}".format(br,vars(br)))  
        
    for t in transformer:
        print("t:{}\n{}".format(t,vars(t)))
        
    for sh in shunt:
        print("sh:{}\n{}".format(sh,vars(sh)))


    start_time = time.time()
    # # # Initialize Solution Vector - V and Q values # # #
    # determine the size of the Y matrix by looking at the total number of nodes in the system
    size_Y = Buses._node_index.__next__()

    # TODO: PART 1, STEP 1 - Complete the function to initialize your solution vector v_init.
    v_init =  np.zeros(size_Y)# create a solution vector filled with zeros of size_Y
    v_init = initialize(size_Y, [*generator,*slack,*load]) # find the initial conditions
    print("initial V", v_init)
    #v_init = np.array([1., 0., 1., 0., 1., 0., 1., 0., 0., 0., 0.])

    # # # Run Power Flow # # #
    powerflow = PowerFlow(case_name, tol, max_iters, enable_limiting, sparse, size_Y)

    # TODO: PART 1, STEP 2 - Complete the PowerFlow class and build your run_powerflow function to solve Equivalent
    #  Circuit Formulation powerflow. The function will return a final solution vector v. Remove run_pf and the if
    #  condition once you've finished building your solver.
    v = powerflow.run_powerflow(v_init, bus, slack, generator, transformer, branch, shunt, load)

    end_time = time.time()
    # # # Process Results # # #
    # TODO: PART 1, STEP 3 - Write a process_results function to compute the relevant results (voltages, powers,
    #  and anything else of interest) and find the voltage profile (maximum and minimum voltages in the case).
    #  You can decide which arguments to pass to this function yourself.
    process_results(v,bus)
    print("total time elapsed: ", end_time - start_time)
    