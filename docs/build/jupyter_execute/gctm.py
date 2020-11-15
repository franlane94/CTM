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

import seaborn as sns
cmap=sns.color_palette('muted')

colours=["black", cmap[4], cmap[1], cmap[6], cmap[0]]

linestyles = ["-", "--", "-.", ":"]

import numpy as np
import matplotlib.pyplot as plt
from gctm import GCTM

# Define the k values

k_vals=np.logspace(-3, 1, 1000)

# Calculate the linear power spectrum at z=0

P_lin_0=GCTM().linear_power(k_vals)

# Calculate the linear power spectrum at z=1

P_lin_1=GCTM().linear_power(k_vals, z_val=1.0)

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear\ at\ }z=0$")
plt.loglog(k_vals, P_lin_1, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Linear\ at\ }z=1$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="upper right", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Calculate the Zel'dovich power spectrum at z=0

P_zel_0=GCTM(nk=1000).zeldovich_power(input_k=k_vals)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_zel_0, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="upper right", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Calculate the Zel'dovich power spectrum at z=0 with kc=5 h/Mpc

P_zel_0_5=GCTM(nk=1000).zeldovich_power(input_k=k_vals, kc=5.0)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_zel_0, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Zel}^\prime\mathrm{dovich}$")
plt.loglog(k_vals, P_zel_0_5, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="lower left", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc

P_gctm_0_5=GCTM(nk=1000).gctm_power(input_k=k_vals, kc=5.0)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_zel_0_5, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Damped\ Zel}^\prime\mathrm{dovich}$")
plt.loglog(k_vals, P_gctm_0_5, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="lower left", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Define redshift values

z_vals=np.linspace(0.0, 200.0, 100)

# Calculate A values

A_vals=np.zeros_like(z_vals)

for i in range(100):

  A_vals[i]=GCTM().linear_growth_factor(z_val=z_vals[i])/GCTM().linear_growth_factor(z_val=99.0)

# Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc with input A

P_gctm_input_A=GCTM(nk=1000).gctm_power(input_k=k_vals, kc=5.0, input_z=z_vals, input_A=A_vals)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_gctm_input_A, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="lower left", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()