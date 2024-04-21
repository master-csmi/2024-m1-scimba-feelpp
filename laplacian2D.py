import feelpp
import feelpp.toolboxes.core as tb

def feelpp_solve_laplacian_2d_toolbox(f=1.0):

    # Create a toolbox environment
    app = tb.Environment(["myapp"], config=tb.localRepository(""))

    # Download the geo file
    geo_path = tb.download( "github:{repo:feelpp,path:feelpp/quickstart/laplacian/cases/feelpp2d/feelpp2d.geo}", worldComm=app.worldCommPtr() )[0]

    # Create a model using the toolbox system
    model = tb.model(dim=2, order=2, filename=geo_path, hsize=0.1)

    # Set the right side function ("f" is the name of the function and str(f) is the expression for the function.)
    model.addFunction("f", str(f))

    # Solve the problem
    model.solve()

    # Export the solution
    model.exportResults()
