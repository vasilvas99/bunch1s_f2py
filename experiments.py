import pybunch1s
import numpy as np
import matplotlib.pyplot as plt
import trajectory_statistics as stats
import integrate_trajectories as itraj

N = 150


def rhs(t, y):
        return pybunch1s.gkrug(y, b=0.5, U=0.05, n=1.0)


y0 = itraj.generate_random_vicinal(N_steps=N, initial_var=0.2)
ts, ys = itraj.integrate_step_trajectory(y0, T_max=600, rhs=rhs)
print(ts[-1])
print("Integration complete")
assert (ts.shape[0] == ys.shape[0])

min_distances = []
for idx, traj in enumerate(ys):
    statistics = stats.l_stat(Y=traj, tk=ts[idx], bdef=1.0)
    min_distances.append(statistics["mind"])

min_distances = np.array(min_distances)
print("Calculating surface statistics (MSI) complete")

l_ts = np.log10(ts[2:])
l_min_distances = np.log10(min_distances[2:])

regr = np.polyfit(l_ts, l_min_distances, 1)
print(regr)

l_y_est = regr[0]*l_ts + regr[1]

plt.xlabel("log10(Time)")
plt.ylabel("log10(Min distance) (MSI)")
plt.title("Min bunch dist (MSI) vs Time for PK")
plt.scatter(l_ts, l_min_distances, label=f"Numerical experiment")
plt.plot(l_ts, l_y_est, color="red", label=f"y = {regr[0]:+.4f} x  {regr[1]:+.4f}")
plt.legend()
plt.show()

