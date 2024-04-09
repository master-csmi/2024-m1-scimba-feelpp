import feelpp as fpp
from feelpp.integrate  import integrate
# import sys


def feelpp_solve_laplacian_2d(f=1.0):
    app = fpp.Environment(["myapp"], config=fpp.localRepository(""))

    geo_path = fpp.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=app.worldCommPtr() )[0]

    mesh = fpp.load(fpp.mesh(dim=2, realdim=2), path=geo_path, size=0.1)

    Xh = fpp.functionSpace(mesh=mesh, family="P2", order=2)
    u = Xh.element()  # Solution function 
    v = Xh.element()  # Test function 

    # Poisson Equation problem
    g = fpp.expr(str(f))  # Constant right side
    a = fpp.form2(Xh)  # Bilinear form
    L = fpp.form1(Xh)  # Linear form
    
    a.integrate(range=fpp.elements(mesh), expr=fpp.gradt(u) * fpp.grad(v))
    L.integrate(range=fpp.elements(mesh), expr=g * v)
    
    # Dirichlet boundary conditions
    a.on(range=fpp.boundaryfaces(mesh), rhs = L, element=u, expr=fpp.expr("0"))
    
    # Solve the Problem
    a.solve(rhs=L, solution=u)
    
    e = fpp.exporter(mesh=mesh, exp=u, name="solution")
    e.add("u",u)
    e.save()

feelpp_solve_laplacian_2d(f=1.0)
