import sys
from scripts.Solve import solve


# path to the grid network RAW file
casename = 'testcases/GS_4_prior_solution.RAW'
casename = 'testcases/IEEE_14_prior_solution.RAW'
casename = 'testcases/IEEE_118_prior_solution.RAW'
#casename = 'testcases/threebus.raw'


# the settings for the solver
settings = {
    "Tolerance": 1E-05,
    "Max Iters": 1000,
    "Sparse": True,
    "Limiting":  True
}

# run the solver
solve(casename, settings)