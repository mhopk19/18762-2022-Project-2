import sys
from scripts.Solve import solve


# path to the grid network RAW file
casename = 'testcases/GS_4_prior_solution.RAW'
#casename = 'testcases/IEEE_14_prior_solution.RAW'
#casename = 'testcases/IEEE_118_prior_solution.RAW'
#casename = 'feasibility-cases/GS_4_stressed.RAW'
#casename = 'testcases/GS_4_prior_solution.RAW'


#casename = 'testcases/PEGASE_9241_flat_start.RAW'
#casename = 'testcases/threebus.raw'


# the settings for the solver
settings = {
    "Tolerance": 1E-07,
    "Max Iters": 100,
    "Sparse": False,
    "Limiting":  False,
    "Feasibility": False
}

# run the solver
solve(casename, settings)