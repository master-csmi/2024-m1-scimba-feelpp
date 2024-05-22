from pathlib import Path

# Import modules
import scimba.equations.domain as domain
import scimba.equations.pde_1d_heat as pde_1d_heat
import scimba.pinns.pinn_tx as pinn_tx
import scimba.pinns.training_tx as training_tx
import scimba.sampling.sampling_ode as sampling_ode
import scimba.sampling.sampling_parameters as sampling_parameters
import scimba.sampling.sampling_pde as sampling_pde
import scimba.sampling.uniform_sampling as uniform_sampling
import torch
import scimba.pinns.pinn_losses as pinn_losses
import scimba.nets.training_tools as training_tools

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch loaded; device is {device}")

torch.set_default_dtype(torch.double)
torch.set_default_device(device)

PI = 3.14159265358979323846


# Test with the Heat equation inside the module PINNs
# where initial and BC are strongly imposed
def main():
    xdomain = domain.SpaceDomain(1,domain.SquareDomain(1, [[0.0, 2.0]]))
    pde = pde_1d_heat.HeatEquation(tdomain=[0.0, 0.03], xdomain=xdomain)
    t_sampler = sampling_ode.TSampler(sampler=uniform_sampling.UniformSampling, ode=pde)
    x_sampler = sampling_pde.XSampler(pde=pde)
    mu_sampler = sampling_parameters.MuSampler(
        sampler=uniform_sampling.UniformSampling, model=pde
    )
    sampler = sampling_pde.PdeTXCartesianSampler(t_sampler, x_sampler, mu_sampler)

    file_name = "test.pth"
    new_training = False #True

    if new_training:
        (
            Path.cwd()
            / Path(training_tx.TrainerPINNSpaceTime.FOLDER_FOR_SAVED_NETWORKS)
            / file_name
        ).unlink(missing_ok=True)

    tlayers = [20, 80, 80, 80, 20, 10]
    network = pinn_tx.MLP_tx(pde=pde, layer_sizes=tlayers, activation_type="sine")
    pinn = pinn_tx.PINNtx(network, pde)

    losses = pinn_losses.PinnLossesData(w_res=1.0)
    optimizers = training_tools.OptimizerData(learning_rate=9e-3,decay=0.99)
    trainer = training_tx.TrainerPINNSpaceTime(
        pde=pde,
        network=pinn,
        losses =losses,
        optimizers=optimizers,
        sampler=sampler,
        file_name=file_name,
        batch_size=4000,
    )

    if new_training:
        trainer.train(epochs=2000, n_collocation=4000)

    trainer.plot(20000,reference_solution=True)


if __name__ == "__main__":
    main()

