import os
import feelpp
from feelpp.toolboxes.cfpdes import *

class Poisson:
  """
  Solves the problem
  -Laplacian u = f   in Omega
  u            = g   in boundary
  
  - with f,g are set by the user
  """
  def __init__(self, dim=2):

    self.dim   = dim
    self.model = dict()
  
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
               plot=None):
    """
    Solves the problem where :
    - h is the mesh size
    - order the polynomial order
    - rhs is the expression of the right-hand side f(x,y)
    """
    self.pb    = cfpdes(dim=self.dim, keyword=f"cfpdes-{self.dim}d")
    self.model = {
      "Name": "Laplacian",
      "ShortName": "Laplacian",
      "Models":
      {
        f"cfpdes-{self.dim}d":
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
        f"cfpdes-{self.dim}d":
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

    feelpp_mesh = feelpp.load(feelpp.mesh(dim=self.dim, realdim=self.dim), fn, h)
    self.pb.setMesh(feelpp_mesh)
    self.pb.setModelProperties(self.model)
    self.pb.init(buildModelAlgebraicFactory=True)
    self.pb.printAndSaveInfo()
    self.pb.solve()
    self.pb.exportResults()
    
    #measures = self.pb.postProcessMeasures().values()
    #return measures
    
    # Plots
    # sudo apt-get update && sudo apt-get install xvfb && pip install --user xvfbwrapper pyvista plotly panel
    if plot != None:
      
      from xvfbwrapper import Xvfb
      import pyvista as pv 

      vdisplay = Xvfb()
      vdisplay.start()
      #pv.set_jupyter_backend('none') 

      mesh = pv.get_reader(f"cfpdes-{self.dim}d.exports/Export.case").read()
      
      pl = pv.Plotter(shape=(1,2))
      pl.add_title(f'Solution P{order}', font_size=18)
      pl.add_mesh(mesh[0], scalars = f"cfpdes.poisson_eq.{name}")
      pl.subplot(0,1)
      pl.add_title('f=' + rhs, font_size=18)
      pl.add_mesh(mesh[0].copy(), scalars = 'cfpdes.expr.rhs')

      pl.link_views()

      if self.dim ==3:
        pl.view_isometric()
      else:
        pl.view_xy()  # if mesh_2D is on the xy plane.        

      pl.screenshot(plot)

    
