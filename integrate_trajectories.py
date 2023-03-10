from numpy import f2py
import numpy as np
import scipy.integrate as sint
import random
import matplotlib.pyplot as plt

N_steps = 50
T_max = 2000
h = 0.1

try:
    import rhs_collection
except Exception as ex:
    with open("rhs-collection.f90") as sourcefile:
        sourcecode = sourcefile.read()
    f2py.compile(sourcecode, modulename='rhs_collection', extension=".f90")
    import rhs_collection

def generate_random_initial_surface(Nsteps):
    list = [i * 1.0 for i in range(1, Nsteps + 1)]
    var = 0.3
    for i in range(Nsteps):
        list[i] = i - 1 + (1.0 - var) + random.uniform(0, 2.0 * var)
    list.sort()
    return np.array(list)

print(rhs_collection.__doc__)

y0 = generate_random_initial_surface(N_steps)
dy = np.zeros_like(y0)

def rhs(t, y):
    return np.squeeze(rhs_collection.g1smm(y, 2.0, 3.0, 2.0))


odeprob = sint.ode(rhs)
odeprob.set_integrator('LSODA', nsteps=3000, rtol=1e-10)
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

plt.show()
