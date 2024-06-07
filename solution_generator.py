import os
import sys
import feelpp
import feelpp.toolboxes.core as tb
from tools.Poisson_feelpp import Poisson_feel

import numpy as np
import pandas as pd
import pyvista as pv

from tools.GmeshRead import mesh2d

def extract_solution(file_path):
    # Fichier .case
    #file_path = '/workspaces/2024-m1-scimba-feelpp/feelppdb/feelpp_cfpde/np_1/cfpdes-2d-p1.exports/Export.case'
    data = pv.read(file_path)

    # Extraire les données de chaque bloc
    for i, block in enumerate(data):
        if block is None:
            continue
        
        print("Champs de points disponibles:", block.point_data.keys())
        print("Champs de cellules disponibles:", block.cell_data.keys())

        solution = block.point_data['cfpdes.poisson.u']
        print("Valeurs de 'cfpdes.poisson.u':")
        print(solution) 


        df = pd.DataFrame(block.point_data)
        print(df.head())
    return solution


def read_nodes(mesh, verbose = True):
    my_mesh = mesh2d(mesh)
    my_mesh.read_mesh()
    coordinates = my_mesh.Nodes
    if verbose :
        print('\nnumber of nodes = ', my_mesh.Nnodes)
        for i in range(my_mesh.Nnodes):
            print(f'\n x_{i} , y_{i} = ', coordinates[i])
    return my_mesh, coordinates

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


# Example usage of the Poisson_feel class


# Define the parameter values you want to iterate over
h_values = [0.1, 0.05]  #, 0.025, 0.0125]  

rhs_values = ['4* pi * (-cos(pi* (x*x + y*y)) + pi * (x*x + y*y)* sin(pi* (x*x + y*y)))', 
              '8*pi*pi*sin(2*pi*x) * sin(2*pi*y)',
              '5/2',
              '-1.0-3*y*x+y*y',
              '4'] 


g_values = ['sin(pi*(x*x + y*y))',
            'sin(2*pi*x) * sin(2*pi*y)', 
            'y',
            '-y*y/2 - x*y*y*y/2 + y*y*y*y/4',
            'x*x/(1+x) + y*y/(1+y)']  


diff_values = ['{1,0,0,1}',
               '{1,0,0,1}',
               '{1,0,0,1}',
               '{1,0,0,1}',
               '{1+x,0,0,1+y}']  

geofile_values = ['omega-2.geo']  

# Create a directory to save the solutions if it doesn't exist
output_dir = 'solutions'
os.makedirs(output_dir, exist_ok=True)

# Iterate over the parameter values
for geofile in geofile_values:
    for h in h_values:
        for rhs, g, diff in zip(rhs_values, g_values, diff_values):
            # Create an instance of P with the current parameters
            P = Poisson_feel(dim=2)
            P(h=h, rhs=rhs, diff=diff, g=g, geofile=geofile)
            
            # Extracting the solution
            solution_path = f"cfpdes-{P.dim}d-p{P.order}.exports/Export.case"
            poisson_u = extract_solution(solution_path)
            
            #  Create a filename for the current solution
            filename = f'solution_h{h}_rhs{rhs.replace("*", "")}_g{g}_geo{os.path.basename(geofile)}.npy'
            filepath = os.path.join(output_dir, filename)
            
            # Save the solution vector
            np.save(filepath, solution_path)

            print(f'Solution saved to {filepath}')

