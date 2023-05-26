import random
import pybunch1s
import numpy as np
import scipy.integrate as sint
import matplotlib.pyplot as plt
from home_integrators import rk23

def generate_random_vicinal(N_steps, initial_var=0.1, bdef=1.0):
    list = [i * 1.0 for i in range(1, N_steps + 1)]
    var = initial_var * bdef
    for i in range(N_steps):
        list[i] = (i - 1) * bdef + (bdef - var) + random.uniform(0, 2.0 * var)
    list.sort()
    return np.array(list)

def integrate_step_trajectory(y0, T_max, rhs, h=0.1, t0=0.0):
    odeprob = sint.ode(rhs)
    odeprob.set_integrator("dopri5", rtol=1e-10)
    odeprob.set_initial_value(y0)

    ts = [t0]
    ys = [y0]
    while ts[-1] <= T_max and odeprob.successful():
        ts.append(odeprob.t + h)
        ys.append(odeprob.integrate(odeprob.t + h))

    ts = np.array(ts)
    ys = np.array(ys)

    return ts, ys


def plot_step_trajectory(ts, ys, N_steps, filepath="initial_results/integration_results.png"):
    plt.rcParams["figure.dpi"] = 400  # set figure dpi

    for i in range(0, N_steps, 1):  # plot every second trajectory
        plt.plot(ts, ys[:, i], linewidth=1)

    plt.savefig(filepath, dpi=500)
    plt.show()

def export_to_txt(ts, ys, filename="step_trajectories.dat"):
    data = np.column_stack((ts, ys))
    header = ["t"] + [f"step_{i}" for i in range(ys.shape[1])]
    header = "\t".join(header)
    np.savetxt(filename, data, delimiter="\t", header=header, comments="")

def main():
    N = 50

    def rhs(t, y):
        # return pybunch1s.g_mm1(y, alpha=0.08, beta=0.0, rho=1.0, n=3.0)
        # return pybunch1s.g_lw(y, p=0, n=2)
        return pybunch1s.g_te(y, p=0, n=2)

    y0 = generate_random_vicinal(N_steps=N, initial_var=0.001, bdef=2)
    # ts, ys = integrate_step_trajectory(y0, T_max=1300, rhs=rhs) # for mm
    # ts, ys = rk23(rhs, y0, T=160000, h0=0.01, tol=0.001)  # for lw
    # export_to_txt(ts, ys)
    plot_step_trajectory(ts, ys, N_steps=N)


if __name__ == "__main__":
    main()
