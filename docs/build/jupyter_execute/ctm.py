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
from ctm import CTM

# Define the k values

k_vals=np.logspace(-3, 1, 1000)

# Calculate the linear power spectrum at z=0

P_lin_0=CTM().linear_power(k_vals)

# Calculate the linear power spectrum at z=1

P_lin_1=CTM().linear_power(k_vals, z_val=1.0)

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear\ at\ }z=0$")
plt.loglog(k_vals, P_lin_1, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Linear\ at\ }z=1$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="upper right", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Calculate the Zel'dovich power spectrum at z=0

P_zel_0=CTM(nk=300).zeldovich_power(input_k=k_vals)

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

P_zel_0_5=CTM(nk=300).zeldovich_power(input_k=k_vals, kc=5.0)

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

P_ctm_0_5=CTM(nk=300).ctm_power(input_k=k_vals, kc=5.0)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_zel_0_5, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Damped\ Zel}^\prime\mathrm{dovich}$")
plt.loglog(k_vals, P_ctm_0_5, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
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

  A_vals[i]=CTM().linear_growth_factor(z_val=z_vals[i])/CTM().linear_growth_factor(z_val=99.0)

# Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc with input A

P_ctm_input_A=CTM(nk=300).ctm_power(input_k=k_vals, kc=5.0, input_z=z_vals, input_A=A_vals)

# Plot the results

plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.loglog(k_vals, P_ctm_input_A, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
plt.legend(loc="lower left", frameon=False, fontsize=14.)
plt.xlim([1e-3, 1])
plt.ylim([1e1, 1e5])
plt.show()

# Compute the linear correlation function
r_lin, corr_lin=CTM().corr_func(k_vals, P_lin_0)

# Compute the Zel'dovich and CTM correlation functions

r_zel, corr_zel=CTM(nk=300).corr_func_zel()
r_ctm, corr_ctm=CTM(nk=300).corr_func_ctm()

# Plot the results

plt.plot(r_lin, r_lin**2*corr_lin, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
plt.plot(r_zel, r_zel**2*corr_zel, color=colours[2], linestyle='--', linewidth=2.2, label=r"$\mathrm{Zel}^\prime\mathrm{dovich}$")
plt.plot(r_ctm, r_ctm**2*corr_ctm, color=colours[4], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Beyond\ Zel}^\prime\mathrm{dovich}$")
plt.xlabel(r"$r\ [\mathrm{Mpc}\ \mathrm{h}^{-1}]$", fontsize=14.)
plt.ylabel(r"$\xi\left(r\right)\ [\mathrm{Mpc}^2\ \mathrm{h}^{-2}]$", fontsize=14.)
plt.legend(loc="upper right", frameon=False, fontsize=14.)
plt.xlim([0, 130])
plt.ylim([-3, 40])
plt.show()