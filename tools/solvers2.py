import os
import feelpp
from feelpp.toolboxes.cfpdes import *
from feelpp.toolboxes.core import *
from feelpp.toolboxes.heat import *
from .lap2D_pinns import Run_laplacian2D, Poisson_2D
from scimba.equations import domain


class HeatSolver:
    def __init__(self, dim=2, order=1):
        self.dim = dim
        self.model = dict()
        self.order = order
        self.pb = None

    def genCube(self, fn, h=0.1):
        """
        Generate a cube geometry following the dimension self.dim
        """
        geo = """SetFactory("OpenCASCADE");
        h={};
        dim={};
        """.format(h, self.dim)

        if self.dim == 2:
            geo += """
            Rectangle(1) = {0, 0, 0, 1, 1, 0};
            Characteristic Length{{ PointsOf{{ Surface{{1}}; }} }} = h;
            Physical Curve("Gamma_D") = {{1,2,3,4}};
            Physical Surface("Omega") = {{1}};
            """
        elif self.dim == 3:
            geo += """
            Box(1) = {0, 0, 0, 1, 1, 1};
            Characteristic Length{{ PointsOf{{ Volume{{1}}; }} }} = h;
            Physical Surface("Gamma_D") = {{1,2,3,4,5,6}};
            Physical Volume("Omega") = {{1}};
            """
        with open(filename, 'w') as f:
            f.write(geo)

    def setModelProperties(self, properties):

        # Create a HeatSolver object
        solver = HeatSolver(dim=3, h=0.1)
        # Define the model properties
        properties = {
        'initial_conditions': ...,  # replace with your initial conditions
        'boundary_conditions': ...,  # replace with your boundary conditions
        'material_properties': ...,  # replace with your material properties
        # ... any other properties ...
        }
        # Set the model properties
        solver.setModelProperties(properties)
        # Now you can call the solve method
        solver.solve()   

    def loadModel(self, model_file):
        with open(model_file, 'r') as f:
            self.model = json.load(f)
            # Add any additional parsing steps here

    def loadMesh(self, mesh_file):
        self.mesh = feelpp.loadMesh(_mesh=feelpp.mesh(dim=self.dim, realdim=self.dim), filename=mesh_file)

    def setupModel(self):
        self.pb = heat(self.mesh)
        self.pb.init()
        self.pb.setModelProperties(self.model)

    def executeSimulation(self, steady=True):
        if steady:
            self.pb.solve()
        else:
            # Implement unsteady simulation here
            # For example, you might use a time-stepping method:
            while not self.pb.timeStepBase().isFinished():
              self.pb.solve()
              self.pb.updateTimeStep()   
            pass   
        
    def postProcess(self):
        self.pb.exportResults()
        # Add any additional post-processing steps here
        # For example, you might calculate and print some statistics about the solution:
        stats = self.pb.statistics()
        print(f"Max temperature: {stats.max()}")
        print(f"Min temperature: {stats.min()}")
        print(f"Average temperature: {stats.mean()}")
        print("Results exported. Open the output file in ParaView or Ensight for visualization.")

        
    def solve(self, solver='feelpp'):
      if solver == 'feelpp':
        print(f"Solving the heat transfer problem...")
            self.pb.setModelProperties(self.model)
            self.pb.init(buildModelAlgebraicFactory=True)
            self.pb.printAndSaveInfo()
            self.pb.solve()
            self.pb.exportResults()
"""
        elif solver == 'scimba':
                print("Solving using Scimba")

                # Create an instance of the problem
                problem = scimba.HeatTransferProblem(self.domain, self.boundary_conditions, self.initial_conditions)

                # Set up the solver
                solver = scimba.Solver(problem)

                # Solve the problem
                solution = solver.solve()

                # Post-process the solution
                self.postProcess(solution)

            # Return the mesh and degrees of freedom
            mesh = pde.mesh()
            function_space = pde.functionSpace()
            dofs = function_space.dofs()
            return mesh, dofs

"""
        print(f"Solving the heat transfer problem for hsize = {self.h}...")
        feelpp_mesh = feelpp.load(feelpp.mesh(dim=self.dim, realdim=self.dim), fn, self.h)
        self.pb = heat(feelpp_mesh)
        self.pb.setModelProperties(self.model)
        self.pb.init(buildModelAlgebraicFactory=True)
        self.pb.printAndSaveInfo()
        self.pb.solve()
        self.pb.exportResults()
        """

# Usage:
solver = HeatSolver(dim=3, order=1)
solver.genCube('my_mesh.geo', h=0.1)
solver.loadModel('my_model.json')
solver.loadMesh('my_mesh.geo')
solver.setupModel()
solver.executeSimulation(steady=False)  # Run an unsteady simulation
solver.postProcess()

class Poisson:
  """
  Solves the problem
  -Laplacian u = f   in Omega
  u            = g   in boundary
  
  - with f,g are set by the user
  """
  def __init__(self, dim=2, order=1):

    self.dim   = dim
    self.model = dict()
    self.order = order
  
  def genCube(self, filename, h=0.1):
    """
    Generate a cube geometry following the dimension  self.dim
    """

    
    geo="""SetFactory("OpenCASCADE");
    h={};
    dim={};
    """.format(h, self.dim)
    
    if self.dim==2 :
        geo+="""
        Rectangle(1) = {0, 0, 0, 1, 1, 0};
        Characteristic Length{ PointsOf{ Surface{1}; } } = h;
        Physical Curve("Gamma_D") = {1,2,3,4};
        Physical Surface("Omega") = {1};
        """
    elif self.dim==3 :
        geo+="""
        Box(1) = {0, 0, 0, 1, 1, 1};
        Characteristic Length{ PointsOf{ Volume{1}; } } = h;
        Physical Surface("Gamma_D") = {1,2,3,4,5,6};
        Physical Volume("Omega") = {1};
        """
    with open(filename, 'w') as f:
        f.write(geo)

  
  
  def __call__(self,
               h,                                       # mesh size 
               order=1,                                 # polynomial order 
               name='Potential',                        # name of the variable u
               rhs='8*pi*pi*sin(2*pi*x)*sin(2*pi*y)',   # right hand side
               diff='{1,0,0,1}',                        # diffusion matrix
               g='0',
               geofile=None,
               plot=None,
               solver='feelpp'):
    """
    Solves the problem where :
    - h is the mesh size
    - order the polynomial order
    - rhs is the expression of the right-hand side f(x,y)
    """
    self.pb    = cfpdes(dim=self.dim, keyword=f"cfpdes-{self.dim}d-p{self.order}")
    self.model = {
      "Name": "Laplacian",
      "ShortName": "Laplacian",
      "Models":
      {
        f"cfpdes-{self.dim}d-p{self.order}":
        {
          "equations":"poisson_eq"
        },
        "poisson_eq":{
          "setup":{
            "unknown":{
              "basis":f"Pch{order}",
              "name":f"{name}",
              "symbol":"u"
            },
            "coefficients":{
              "c": f"{diff}:x:y" if self.dim == 2 else f"{diff}:x:y:z",
              "f": f"{rhs}:x:y"  if self.dim == 2 else f"{rhs}:x:y:z"
            }
          }
        }
      },
      "Materials":
      {
        "Omega":
        {
          "markers":["Omega"]
        }
      },
      "BoundaryConditions":
      {
        "poisson_eq":
        {
          "Dirichlet":
          {
            "g":
            {
              "markers":["Gamma_D"],
              "expr":f"{g}:x:y"
            }
          }
        }
      },
      "PostProcess":
      {
        f"cfpdes-{self.dim}d-p{self.order}":
        {
          "Exports":
          {
            "fields":["all"],
            "expr":{
              "rhs"    : f"{rhs}:x:y" if self.dim == 2 else f"{rhs}:x:y:z"
            }
          }
        }
      }
    }

    fn = None
    if geofile is None:
      fn = f'omega-{self.dim}.geo'
      self.genCube(fn, h)
    else:
      fn = geofile
    
        
##________________________

    if solver == 'feelpp':
      print(f"Solving the laplacian problem for hsize = {h}...")
      feelpp_mesh = feelpp.load(feelpp.mesh(dim=self.dim, realdim=self.dim), fn, h)
      self.pb.setMesh(feelpp_mesh)
      self.pb.setModelProperties(self.model)
      self.pb.init(buildModelAlgebraicFactory=True)
      self.pb.printAndSaveInfo()
      self.pb.solve()
      self.pb.exportResults()
      #measures = self.pb.postProcessMeasures().values()
      """
      try:
        import pandas as pd
        df=pd.DataFrame([measures])
        print(df)
      except ImportError:
        print("cannot import pandas, no problem it was just a test")
    """
##______________________
      #mesh = self.pb.mesh()
      #function_space = self.pb.functionSpace()
      #dofs = function_space.dofs()
##________________________
      #return mesh, dofs
##________________________

    elif solver == 'scimba':
      # Placeholder for Scimba solving logic
      print("Solving using Scimba")

      # Define a square domain
      xdomain = domain.SpaceDomain(2, domain.SquareDomain(2, [[0.0, 1.0], [0.0, 1.0]]))

      # Create an instance of the Poisson problem
      pde = Poisson_2D(xdomain)

      # Run the training
      Run_laplacian2D(pde)
##________________________

    # Plots
    # sudo apt-get update && sudo apt-get install xvfb && pip install --user xvfbwrapper pyvista plotly panel
    if plot != None:
      
      from xvfbwrapper import Xvfb
      import pyvista as pv 

      #vdisplay = Xvfb()
      #vdisplay.start()
      pv.set_jupyter_backend('static') 
      pv.start_xvfb()

      mesh = pv.get_reader(f"cfpdes-{self.dim}d-p{self.order}.exports/Export.case").read()
      
      pl = pv.Plotter(shape=(1,2))
      pl.add_title(f'Solution P{order}', font_size=18)
      pl.add_mesh(mesh[0], scalars = f"cfpdes.poisson_eq.{name}")
      pl.subplot(0,1)
      pl.add_title('f=' + rhs, font_size=18)
      pl.add_mesh(mesh[0].copy(), scalars = 'cfpdes.expr.rhs')
      pl.show()
      pl.link_views()

      if self.dim ==3:
        pl.view_isometric()
      else:
        pl.view_xy()  # if mesh_2D is on the xy plane.        
      
      pl.show()
      pl.screenshot(plot)

    
