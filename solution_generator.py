import os
import sys
import feelpp
import feelpp.toolboxes.core as tb
from tools.Poisson_feelpp import Poisson_feel

# mandatory things
sys.argv = ["feelpp_app"]
e = feelpp.Environment(sys.argv,
                       opts=tb.toolboxes_options("coefficient-form-pdes", "cfpdes"),
                       config=feelpp.localRepository('feelpp_cfpde'))

# ------------------------------------------------------------------------- #
# Poisson problem
# - div (diff * grad (u)) = f    in Omega
#                     u   = g    in Gamma_D
# Omega = domain, either cube or ball
# Approx = lagrange Pk of order order
# mesh of size h


# Example usage of the Poisson class
# Create an instance of the Poisson class for a 2-dimensional problem
P = Poisson_feel(dim=2)

# Solve the Poisson problem with the specified parameters
P(h=0.08, rhs='-1.0-1*y*x+y*y', g='0', geofile='geo/disk.geo')

import numpy as np
import pandas as pd

# Define the parameter values you want to iterate over
h_values = [0.08]  # You can add more values if needed
rhs_values = ['-1.0-1*y*x+y*y']  # You can add more expressions if needed
g_values = ['0']  # You can add more expressions if needed
geofile_values = ['geo/disk.geo']  # You can add more file paths if needed

# Create a directory to save the solutions if it doesn't exist
output_dir = 'solutions'
os.makedirs(output_dir, exist_ok=True)

# Iterate over the parameter values
for h in h_values:
    for rhs in rhs_values:
        for g in g_values:
            for geofile in geofile_values:
                # Create an instance of P with the current parameters
                p = P(h=h, rhs=rhs, g=g, geofile=geofile)
                
                # Solve the problem and get the solution vector
                solution_vector = p.solve()
                
                # Create a filename for the current solution
                filename = f'solution_h{h}_rhs{rhs.replace("*", "")}_g{g}_geo{os.path.basename(geofile)}.csv'
                filepath = os.path.join(output_dir, filename)
                
                # Save the solution vector to a CSV file
                np.savetxt(filepath, solution_vector, delimiter=',')

                print(f'Solution saved to {filepath}')

