import pybunch1s
import numpy as np
import matplotlib.pyplot as plt
import trajectory_statistics as stats
import integrate_trajectories as itraj

N = 150


def rhs(t, y):
    return pybunch1s.g_mm1(y, alpha=0.08, beta=0.0, rho=1.0, n=1.0)


y0 = itraj.generate_random_vicinal(N_steps=N, initial_var=0.2)
ts, ys = itraj.integrate_step_trajectory(y0, T_max=1300, rhs=rhs)

print("Integration complete")
assert (ts.shape[0] == ys.shape[0])

min_distances = []
for idx, traj in enumerate(ys):
    statistics = stats.l_stat(Y=traj, tk=ts[idx], bdef=1.0)
    min_distances.append(statistics["mind"])

min_distances = np.array(min_distances)
print("Calculating surface statistics (MSI) complete")

l_ts = np.log(ts[2:])
l_min_distancces = np.log(min_distances[2:])

regr = np.polyfit(l_ts, l_min_distancces, 1)
print(regr)


plt.xlabel("Log(Time)")
plt.ylabel("Log(Min distance) (MSI)")
plt.title("Min bunch dist (MSI) vs Time for gmm1")
plt.plot(l_ts, l_min_distancces)
plt.show()

