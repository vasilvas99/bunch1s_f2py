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

    #   write(10,3075) h+1,w,w/h,mind,xk

p = Path(sys.argv[1]).resolve(strict=True)

try:
    second_arg = sys.argv[2]
    if second_arg == "--reverse-cols":
        col_names = ["N", "w", "w/h", "mind", "time"]
    else:
        print(f"Second CLI arg does not make sense. Expected: --reverse-cols, got: {second_arg}. Using default column ordering")
        col_names = ["time", "N", "w", "w/h", "mind"]
except IndexError as _:
        col_names = ["time", "N", "w", "w/h", "mind"]

print(f"Using column ordering for MSI output: {col_names=}")

data = pd.read_fwf(p, header=None ,names=col_names)
data_log = np.log10(data)

data_log.replace([np.inf, -np.inf], np.nan, inplace=True)
data_log.dropna(inplace=True)


PLOT_KEY = "N"
time_log = data_log["time"]
bunch_statistic_log = data_log[PLOT_KEY]

print(time_log)

slope, intercept = np.polyfit(time_log, bunch_statistic_log, deg=1)
bunch_statistic_log_fit = slope*time_log + intercept

plt.scatter(time_log, bunch_statistic_log, color="red", label="MSI data")
plt.plot(time_log, bunch_statistic_log_fit, label = line_str(slope, intercept))

plt.xlim([np.min(time_log), np.max(time_log)])
plt.ylim([np.min(bunch_statistic_log), np.max(bunch_statistic_log)])

plt.xlabel(r"$\log(t)$")
plt.ylabel(f"$\\log({PLOT_KEY})$")
plt.legend()

plt.show()