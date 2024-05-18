# This has to be executed at every time we run the docker image (could be automatized ?)
# sudo apt-get update && sudo apt-get install xvfb && pip install --user xvfbwrapper pyvista plotly panel
# ./tools/load_xvfb.sh  for plotting in docker term


import sys
import feelpp
import feelpp.toolboxes.core as tb

from tools.solvers import Poisson

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

# 2D with varying RHS
P = Poisson(dim = 2)
P(h=0.08,  rhs='-1.0-1*y*x+y*y', g='0', order=1, geofile='geo/disk.geo', plot='2d.png')
P(h=0.08,  rhs='-1.0-1*y*x+y*y', g='0', order=1, geofile='geo/disk.geo', solver='scimba')


P(h=0.1,  rhs='-1.0-2*y*x+y*y', g='0', order=1, plot='f2.png')
P(h=0.1,  rhs='-1.0-2*y*x+y*y', g='0', order=1, solver ='scimba')

P(h=0.1,  rhs='-1.0-3*y*x+y*y', g='y', order=1, plot='f3.png')
P(h=0.1,  rhs='-1.0-3*y*x+y*y', g='y', order=1, solver ='scimba')


P(h=0.1,  rhs='-1.0-4*y*x+y*y', g='x', order=1, plot='f4.png')
P(h=0.1,  rhs='-1.0-4*y*x+y*y', g='x', order=1, solver ='scimba')

P(h=0.05, rhs='1',              g='0', order=1, plot='f5.png')
P(h=0.05, rhs='1 + x-x',              g='0', order=1, solver ='scimba')

# # 2D with varying anisotropy
#P = Poisson(dim = 2)
P(h=0.1, diff='{1.0,0,0,x*y}', rhs='1', plot='d1.png')
P(h=0.1, diff='(1+x,0,0,1+y)', rhs='1+ x-x', solver='scimba')

P(h=0.1, diff='{x,y,-y,x+y}',  rhs='1', plot='d3.png')
P(h=0.1, diff='(x,y,-y,x+y)',  rhs='1+ x-x', solver='scimba')

# 3D with varying anisotropy and non-homogeneous Dirichlet BC
#P = Poisson(dim = 3)
#P(h=0.08, diff='{1,0,0,0,x+1,0,0,0,1+x*y}', g = 'x', rhs='x*y*z', geofile = 'geo/cube.geo', plot='3d.png') # au travers __call__



# ce qu'on voudrait

# 1 pouvoir résoudre dans scimba juste en changeant le flag "solver"
#P(h=0.08, diff='{1,0,0,0,x+1,0,0,0,1+x*y}', g = 'x', rhs='x*y*z', geofile = 'geo/cube.geo', plot='3d.png', solver='scimba') # au travers __call__

# pouvoir rcupérer le mesh et les dofs en sortie de P(...)
#mean


"""
pour scimba
# eval('__') rhs = x1 
f = eval(rhs)
rhs.replace('x','x1')

"""
#mesh, dofs = P(...)....
#utiliser diff sur scimba dx*()+dy(+dz) cas matrice cte