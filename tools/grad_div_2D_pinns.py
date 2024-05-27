from pathlib import Path

import scimba.nets.training_tools as training_tools
import scimba.pinns.pinn_losses as pinn_losses
import scimba.pinns.pinn_x as pinn_x
import scimba.pinns.training_x as training_x
import scimba.sampling.sampling_parameters as sampling_parameters
import scimba.sampling.sampling_pde as sampling_pde
import scimba.sampling.uniform_sampling as uniform_sampling
import torch
from scimba.equations import domain
from scimba.equations.pdes import AbstractPDEx

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch loaded; device is {device}")

torch.set_default_dtype(torch.double)
torch.set_default_device(device)

PI = 3.14159265358979323846


class GradDiv(AbstractPDEx):
    r"""

    .. math::

        \nabla (\nabla \cdot u) + u = f

    """

    def __init__(self, space_domain):
        super().__init__(
            nb_unknowns=2,
            space_domain=space_domain,
            nb_parameters=1,
            parameter_domain=[[0.75, 0.75 + 1e-4]],
        )

        self.first_derivative = True
        self.second_derivative = True

    def make_data(self, n_data):
        pass

    def bc_residual(self, w, x, mu, **kwargs):
        u1, u2 = self.get_variables(w)
        return u1, u2

    def residual(self, w, x, mu, **kwargs):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        u1, u2 = self.get_variables(w)
        u1_xx, u2_xx = self.get_variables(w, "w_xx")
        u1_xy, u2_xy = self.get_variables(w, "w_xy")
        u1_yy, u2_yy = self.get_variables(w, "w_yy")

        cos_x1 = torch.cos(2.0 * PI * x1)
        cos_x2 = torch.cos(2.0 * PI * x2)
        sin_x1 = torch.sin(2.0 * PI * x1)
        sin_x2 = torch.sin(2.0 * PI * x2)

        f1 = (1 - 4 * PI**2) * sin_x1 * sin_x2 + 4 * PI**2 * alpha * cos_x1 * cos_x2
        f2 = (1 - 4 * PI**2 * alpha) * sin_x1 * sin_x2 + 4 * PI**2 * cos_x1 * cos_x2

        return u1_xx + u2_xy + u1 - f1, u1_xy + u2_yy + u2 - f2

    def post_processing(self, x, mu, w):
        x1, x2 = x.get_coordinates()
        return x1 * (1 - x1) * x2 * (1 - x2) * w

    def reference_solution(self, x, mu):
        x1, x2 = x.get_coordinates()
        alpha = self.get_parameters(mu)
        return torch.cat(
            (
                torch.sin(2.0 * PI * x1) * torch.sin(2.0 * PI * x2),
                alpha * torch.sin(2.0 * PI * x1) * torch.sin(2.0 * PI * x2),
            ),
            axis=1,
        )


def main():
    xdomain = domain.SpaceDomain(2, domain.SquareDomain(2, [[0.0, 1.0], [0.0, 1.0]]))
    pde = GradDiv(xdomain)
    x_sampler = sampling_pde.XSampler(pde=pde)
    mu_sampler = sampling_parameters.MuSampler(
        sampler=uniform_sampling.UniformSampling, model=pde
    )
    sampler = sampling_pde.PdeXCartesianSampler(x_sampler, mu_sampler)

    file_name = "test.pth"
    new_training = True
    # new_training = False

    train = True

    if new_training:
        (
            Path.cwd()
            / Path(training_x.TrainerPINNSpace.FOLDER_FOR_SAVED_NETWORKS)
            / file_name
        ).unlink(missing_ok=True)

    tlayers = [20, 20, 20, 20, 20]
    network = pinn_x.DisMLP_x(pde=pde, layer_sizes=tlayers, activation_type="sine")
    pinn = pinn_x.PINNx(network, pde)

    losses = pinn_losses.PinnLossesData(w_res=1.0)
    optimizers = training_tools.OptimizerData(
        learning_rate=1.5e-2, decay=0.99, switch_to_LBFGS=True
    )
    trainer = training_x.TrainerPINNSpace(
        pde=pde,
        network=pinn,
        losses=losses,
        optimizers=optimizers,
        sampler=sampler,
        file_name=file_name,
        batch_size=5000,
    )

    if train:
        trainer.train(epochs=2000, n_collocation=5000, n_bc_collocation=0, n_data=0)

    trainer.plot(20000, reference_solution=True)


if __name__ == "__main__":
    main()
