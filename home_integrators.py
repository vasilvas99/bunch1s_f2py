import numpy as np
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

def rk23(rhs, y0, T, h0, tol, t0=0):
    # local truncation error estimator
    rk23_te = lambda k1, k2, k3: np.linalg.norm((1 / 3) * (-k1 - k2 + 2 * k3), np.inf)
    y0 = np.atleast_1d(y0)

    t = [t0]
    y = [y0]
    h = h0
    pbar = tqdm(desc="RK23 Integration", total = T, unit = "simul step")
    while t[-1] <= T:
        k1 = h * rhs(t[-1], y[-1])
        k2 = h * rhs(t[-1] + h, y[-1] + k1)
        k3 = h * rhs(t[-1] + h / 2, y[-1] + k1 / 4 + k2 / 4)

        while rk23_te(k1, k2, k3) > tol:
            h *= (tol / rk23_te(k1, k2, k3)) ** (1 / 3)
            k1 = h * rhs(t[-1], y[-1])
            k2 = h * rhs(t[-1] + h, y[-1] + k1)
            k3 = h * rhs(t[-1] + h / 2, y[-1] + k1 / 4 + k2 / 4)

        y_new = y[-1] + k1 / 6 + k2 / 6 + 2 * k3 / 3
        t_new = t[-1] + h

        y.append(y_new)
        t.append(t_new)
        pbar.update(h)
        h *= (tol / rk23_te(k1, k2, k3)) ** (1 / 3)
    pbar.close()
    return np.array(t), np.squeeze(np.array(y))