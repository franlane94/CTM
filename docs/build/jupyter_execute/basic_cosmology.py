import matplotlib as mpl

mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
mpl.rcParams["text.usetex"] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cm'
mpl.rcParams["lines.linewidth"] = 2.2
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.labelsize"] = 14.
mpl.rcParams["xtick.top"] = True
mpl.rcParams["xtick.labelsize"] = 14.
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.right"] = True
mpl.rcParams["ytick.labelsize"] = 14.
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.minor.bottom"] = False
mpl.rcParams["xtick.minor.top"] = False
mpl.rcParams["ytick.minor.left"] = False
mpl.rcParams["ytick.minor.right"] = False

linestyles = ["-", "--", "-.", ":"]

import numpy as np
import matplotlib.pyplot as plt
from ctm import Cosmo

# Define the redshift values

z_vals=np.linspace(0.0, 200.0, 100)

# Calculate the linear growth factor values

D_1_vals=Cosmo(h=0.7, omega0_b=0.02233, omega0_cdm=0.112, n_s=0.96, sigma_8=0.8, k_max=10.0, verbose=False, gauge='sync', output='mPk').calc_linear_growth(z_vals)

# Plot the results

plt.plot(z_vals, D_1_vals, color="black", linestyle='-', linewidth=2.2, alpha=0.8)
plt.xlabel(r"$z$", fontsize=14.)
plt.ylabel(r"$D_1$", fontsize=14.)
plt.xlim([-2, 202])
plt.ylim([-0.2, 1.2])
plt.show()

# Calculate the independent growth factor values

f_vals=Cosmo(h=0.7, omega0_b=0.02233, omega0_cdm=0.112, n_s=0.96, sigma_8=0.8, k_max=10.0, verbose=False, gauge='sync', output='mPk').calc_independent_linear_growth(z_vals)

# Plot the results

plt.plot(z_vals, f_vals, color="black", linestyle='-', linewidth=2.2, alpha=0.8)
plt.xlabel(r"$z$", fontsize=14.)
plt.ylabel(r"$f$", fontsize=14.)
plt.xlim([-2, 202])
plt.ylim([0.5, 1.2])
plt.show()