from pathlib import Path

import scimba.nets.training_tools as training_tools
import scimba.pinns.pinn_losses as pinn_losses
import scimba.pinns.pinn_x as pinn_x
import scimba.pinns.training_x as training_x
import scimba.sampling.sampling_parameters as sampling_parameters
import scimba.sampling.sampling_pde as sampling_pde
import scimba.sampling.uniform_sampling as uniform_sampling
import torch
from scimba.equations import domain, pdes


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch loaded; device is {device}")

torch.set_default_dtype(torch.double)
torch.set_default_device(device)

PI = 3.14159265358979323846
ELLIPSOID_A = 4 / 3
ELLIPSOID_B = 1 / ELLIPSOID_A


class PoissonDisk2D(pdes.AbstractPDEx):
    def __init__(self, space_domain, 
                 rhs='4*pi*sin(pi*(x*x + y*y)) - 4*pi*pi*(x*x + y*y)*cos(pi*(x*x + y*y))', 
                 diff='(1,0,0,1)', 
                 g='0',
                 u_exact = 'sin(pi*(x*x + y*y))'):
        super().__init__(
            nb_unknowns=1,
            space_domain=space_domain,
            nb_parameters=1,
            parameter_domain=[[0.5, 1]],
        )

        self.rhs = rhs
        self.diff = diff
        self.g = g
        self.u_exact = u_exact 
        self.first_derivative = True
        self.second_derivative = True

    def bc_residual(self, w, x, mu, **kwargs):
        u = self.get_variables(w)
        x1, x2 = x.get_coordinates()
        g_evaluated = eval(self.g, {'x': x1, 'y': x2, 'pi': PI, 'sin' : torch.sin, 'cos': torch.cos})
        return u - g_evaluated

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        u_xx = self.get_variables(w, "w_xx")
        u_yy = self.get_variables(w, "w_yy")

        f = eval(self.rhs, {'x': x1, 'y': x2, 'pi': PI, 'sin': torch.sin, 'cos': torch.cos})
        diff = eval(self.diff, {'x': x1, 'y': x2, 'pi': PI, 'sin' : torch.sin, 'cos': torch.cos})
        
        return u_xx * diff[0] + u_yy * diff[3] + f  


    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        return eval(self.u_exact, {'x': x1, 'y': x2, 'pi': PI, 'sin': torch.sin, 'cos': torch.cos})


class Poisson_2D(pdes.AbstractPDEx):
    def __init__(self, space_domain,  
                 rhs='8*pi*pi*sin(2*pi*x)*sin(2*pi*y)', 
                 diff='(1,0,0,1)', 
                 g='0',
                 u_exact='sin(2 * pi * x) * sin(2 * pi * y)'):
        super().__init__(
            nb_unknowns=1,
            space_domain=space_domain,
            nb_parameters=1,
            parameter_domain=[[0.50000, 0.500001]],           
        )
        self.rhs = rhs
        self.diff = diff
        self.g = g
        self.u_exact = u_exact
        self.first_derivative = True
        self.second_derivative = True

    def bc_residual(self, w, x, mu, **kwargs):
        u = self.get_variables(w)
        # Évaluation de la condition aux limites g
        x1, x2 = x.get_coordinates()
        g_evaluated = eval(self.g, {'x': x1, 'y': x2, 'pi': PI, 'sin' : torch.sin, 'cos': torch.cos})
        return u - g_evaluated

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        u_xx = self.get_variables(w, "w_xx")
        u_yy = self.get_variables(w, "w_yy")

        f = eval(self.rhs, {'x': x1, 'y': x2, 'pi': PI, 'sin' : torch.sin, 'cos': torch.cos})
        diff = eval(self.diff, {'x': x1, 'y': x2, 'pi': PI, 'sin' : torch.sin, 'cos': torch.cos})
    
        return u_xx * diff[0] + u_yy * diff[3] + f  

    """
    def post_processing(self, x, mu, w):
        x1, x2 = x.get_coordinates()
        return x1 * (1 - x1) * x2 * (1 - x2) * w
    """
    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        return eval(self.u_exact, {'x': x1, 'y': x2, 'pi': PI, 'sin': torch.sin, 'cos': torch.cos})

class Poisson_2D_ellipse(pdes.AbstractPDEx):
    def __init__(self, space_domain):
        super().__init__(
            nb_unknowns=1,
            space_domain=space_domain,
            nb_parameters=0,
            parameter_domain=[[0.99999, 1]],
        )

        self.first_derivative = True
        self.second_derivative = True

    def make_data(self, n_data):
        pass

    def bc_residual(self, w, x, mu, **kwargs):
        u = self.get_variables(w)
        return u

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        u_xx = self.get_variables(w, "w_xx")
        u_yy = self.get_variables(w, "w_yy")
        f = 1
        return u_xx + u_yy + f

    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        x1_0, x2_0 = self.space_domain.large_domain.center
        a, b = ELLIPSOID_A, ELLIPSOID_B
        rho = 0.5 / (1 / a**2 + 1 / b**2)
        return rho * (1 - ((x1 - x1_0) / a) ** 2 - ((x2 - x2_0) / b) ** 2)


def disk_to_ellipse(x):
    x1, x2 = (x[:, i, None] for i in range(2))
    return torch.cat((x1 * ELLIPSOID_A, x2 * ELLIPSOID_B), axis=1)


def Jacobian_disk_to_ellipse(x):
    x1, x2 = (x[:, i, None] for i in range(2))
    return ELLIPSOID_A, 0, 0, ELLIPSOID_B


def disk_to_potato(x):
    x1, x2 = (x[:, i, None] for i in range(2))
    x = x1 - 0.5 * x2**2 + 0.3 * torch.sin(x2)
    y = x2 + 0.1 * x + 0.12 * torch.cos(x)
    return torch.cat((x, y), axis=1)


def Jacobian_disk_to_potato(x):
    x1, x2 = (x[:, i, None] for i in range(2))
    raise ValueError("Jacobian_disk_to_potato is not implemented")
    return 0, 0, 0, 0


def Run_laplacian2D(pde, bc_loss_bool=False, w_bc=0, w_res=1.0):
    x_sampler = sampling_pde.XSampler(pde=pde)
    mu_sampler = sampling_parameters.MuSampler(
        sampler=uniform_sampling.UniformSampling, model=pde
    )
    sampler = sampling_pde.PdeXCartesianSampler(x_sampler, mu_sampler)

    file_name = "test.pth"
    #new_training = False
    new_training = True

    if new_training:
        (
            Path.cwd()
            / Path(training_x.TrainerPINNSpace.FOLDER_FOR_SAVED_NETWORKS)
            / file_name
        ).unlink(missing_ok=True)

    tlayers = [20, 20, 20, 20, 20]
    network = pinn_x.MLP_x(pde=pde, layer_sizes=tlayers, activation_type="sine")
    pinn = pinn_x.PINNx(network, pde)
    losses = pinn_losses.PinnLossesData(
        bc_loss_bool=bc_loss_bool, w_res=w_res, w_bc=w_bc
    )
    optimizers = training_tools.OptimizerData(learning_rate=1.2e-2, decay=0.99)
    trainer = training_x.TrainerPINNSpace(
        pde=pde,
        network=pinn,
        sampler=sampler,
        losses=losses,
        optimizers=optimizers,
        file_name=file_name,
        batch_size=5000,
    )

    if not bc_loss_bool:
        if new_training:
            trainer.train(epochs=100, n_collocation=5000, n_data=0)
    else:
        if new_training:
            trainer.train(
                epochs=100, n_collocation=5000, n_bc_collocation=1000, n_data=0
            )

    trainer.plot(20000, reference_solution=True)
    # trainer.plot_derivative_mu(n_visu=20000)
    return network, pde

def solution_array(pde):
    network, pde = Run_laplacian2D(pde)
    # Extract solution function u
    u = network.forward


if __name__ == "__main__":
    # Laplacien strong Bc on Square with nn
    xdomain = domain.SpaceDomain(2, domain.SquareDomain(2, [[0.0, 1.0], [0.0, 1.0]]))
    print(xdomain)

    u_exact='sin(2 * pi * x) * sin(2 * pi * y)'       
    pde = Poisson_2D(xdomain,  rhs='8*pi*pi*sin(2*pi*x)*sin(2*pi*y)', g='0',  u_exact=u_exact)
    network, pde = Run_laplacian2D(pde)

    """
    pde = Poisson_2D(xdomain, rhs='-1.0-4*y*x+y*y', g='x')
    network, pde = Run_laplacian2D(pde)
    
    xdomain = domain.SpaceDomain(2, domain.DiskBasedDomain(2, center=[0.0, 0.0], radius=1.0))
    u_exact =  'sin(pi*(x*x + y*y))'
    rhs = '4*pi*sin(pi*(x*x + y*y)) - 4*pi*pi*(x*x + y*y)*cos(pi*(x*x + y*y))'

    pde_disk = PoissonDisk2D(xdomain,  rhs= rhs, g= '0', u_exact=u_exact)
    network, pde = Run_laplacian2D(pde_disk)
    """

    # Extract solution function u
    u = network.forward
    # Generate example input data and get the solution from the network
    input_tensor = torch.randn(100, 2)  # Generate some example input data
    mu = torch.tensor([[0.5]]).repeat(input_tensor.size(0), 1)  # Expand and repeat mu to match the batch size
    solution_tensor = u(input_tensor, mu)   # Get the solution from the network

    # Convert the tensor to a NumPy array
    solution_array = solution_tensor.detach().numpy()
    print('solution = ', solution_array)