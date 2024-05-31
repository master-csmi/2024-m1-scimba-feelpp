# common.py
import torch
from scimba.equations import domain, pdes

def Run_laplacian2D(pde, bc_loss_bool=False, w_bc=0, w_res=1.0):
    import scimba.nets.training_tools as training_tools
    import scimba.pinns.pinn_losses as pinn_losses
    import scimba.pinns.pinn_x as pinn_x
    import scimba.pinns.training_x as training_x
    import scimba.sampling.sampling_parameters as sampling_parameters
    import scimba.sampling.sampling_pde as sampling_pde
    import scimba.sampling.uniform_sampling as uniform_sampling
    from pathlib import Path

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.set_default_dtype(torch.double)
    torch.set_default_device(device)

    x_sampler = sampling_pde.XSampler(pde=pde)
    mu_sampler = sampling_parameters.MuSampler(
        sampler=uniform_sampling.UniformSampling, model=pde
    )
    sampler = sampling_pde.PdeXCartesianSampler(x_sampler, mu_sampler)

    file_name = "test.pth"
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
            trainer.train(epochs=600, n_collocation=5000, n_data=0)
    else:
        if new_training:
            trainer.train(
                epochs=600, n_collocation=5000, n_bc_collocation=1000, n_data=0
            )

    trainer.plot(20000, reference_solution=True)
    return network, pde

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
        g_evaluated = eval(self.g, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin' : torch.sin, 'cos': torch.cos})
        return u - g_evaluated

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        u_xx = self.get_variables(w, "w_xx")
        u_yy = self.get_variables(w, "w_yy")

        f = eval(self.rhs, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin': torch.sin, 'cos': torch.cos})
        diff = eval(self.diff, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin' : torch.sin, 'cos': torch.cos})
        
        return u_xx * diff[0] + u_yy * diff[3] + f  

    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        return eval(self.u_exact, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin': torch.sin, 'cos': torch.cos})

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
        g_evaluated = eval(self.g, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin' : torch.sin, 'cos': torch.cos})
        return u - g_evaluated

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        u_xx = self.get_variables(w, "w_xx")
        u_yy = self.get_variables(w, "w_yy")

        f = eval(self.rhs, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin' : torch.sin, 'cos': torch.cos})
        diff = eval(self.diff, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin' : torch.sin, 'cos': torch.cos})
    
        return u_xx * diff[0] + u_yy * diff[3] + f  

    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        return eval(self.u_exact, {'x': x1, 'y': x2, 'pi': 3.14159265358979323846, 'sin': torch.sin, 'cos': torch.cos})
