import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

def line_str(slope, intercept):
    if intercept > 0:
        return f"y = {slope:.4f} x  + {intercept:.4f}"
    else:
        return f"y = {slope:.4f} x  - {np.abs(intercept):.4f}"


p = Path(sys.argv[1]).resolve(strict=True)
data = pd.read_fwf(p, header=None ,names=["time", "h + 1", "w", "w/h", "mind"])

time = data["time"]
bunch_width = data["w"]

mask = time > 0
time = time[mask]
bunch_width = bunch_width[mask]

t_log = np.log(time)
w_log = np.log(bunch_width)

slope, intercept = np.polyfit(t_log, w_log, deg=1)


w_log_fit = slope*t_log + intercept

plt.scatter(t_log, w_log, color="red", label="MSI data")
plt.xlim([np.min(t_log), np.max(t_log)])
plt.ylim([np.min(w_log), np.max(w_log)])
plt.xlabel(r"$\log(t)$")
plt.ylabel(r"$\log(w)$")
plt.plot(t_log, w_log_fit, label = line_str(slope, intercept))
plt.legend()
plt.show()