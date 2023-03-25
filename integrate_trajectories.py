import numpy as np
import scipy.integrate as sint
import random
import matplotlib.pyplot as plt
import pybunch1s

N_steps = 50
T_max = 2000
h = 0.1

def generate_random_initial_surface(Nsteps):
    list = [i * 1.0 for i in range(1, Nsteps + 1)]
    var = 0.3
    for i in range(Nsteps):
        list[i] = i - 1 + (1.0 - var) + random.uniform(0, 2.0 * var)
    list.sort()
    return np.array(list)

y0 = generate_random_initial_surface(N_steps)
dy = np.zeros_like(y0)

def rhs(t, y):
    return np.squeeze(pybunch1s.g1smm(y, 2.0, 3.0, 2.0))


odeprob = sint.ode(rhs)
odeprob.set_integrator('LSODA', rtol=1e-10)
odeprob.set_initial_value(y0)


ts = [0]
ys = [y0]
while odeprob.successful() and odeprob.t < T_max:
    ts.append(odeprob.t+h)
    ys.append(odeprob.integrate(odeprob.t+h))

ts = np.array(ts)
ys = np.array(ys)

plt.rcParams['figure.dpi'] = 400 # set figure dpi

for i in range(0, N_steps, 2): # plot every second trajectory
    plt.plot(ts, ys[:,i], linewidth=1)
    
plt.savefig("initial_results/g1smm.png", dpi=300)
plt.show()
