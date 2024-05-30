import scimba
import feelpp

import toolbox_rh

# Create a uniform mesh  of100 points in R^2
m = Mesh(100) # ....

# sketch of Solving  with Scimba
u_scimba = scimba.pinn_x.poisson(f='...', m)

# sketch of solving with Feel++ CFPDES
u_feelpp = toolbox_rh.solve_Poisson(f='sin(x)*sin(y)', m)

#
err = norm(u_scimba - u_feelpp)



