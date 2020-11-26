Using the GCTM in real space
=============================

Documention introducing the GCTM class

.. py:class:: GCTM(min_k=1e-7, max_k=2e3, min_k_calc=1e-5, max_k_calc=100.0, nk=3000, h=0.6737, omega0_b=0.02233, omega0_cdm=0.11977, n_s=0.9652, sigma_8=0.8101, verbose=False, gauge='sync', output='mPk', **kwargs)

Class to calculate the linear, Zel'dovich and GCTM power spectra in real and redshift space

:Parameters: - min_k (*int*) : the minimum :math:`k` value used to calculate loop integrals (**default:       min_k=1e-7**)
             - max_k (*int*) : the maximum :math:`k` value used to calculate loop integrals (**default: max_k=2e3**)
             - min_k_calc (*int*) : the minimum :math:`k` value used to calculate the power spectrum (**default: min_k_calc=1e-5**)
             - max_k_calc (*int*) : the maximum :math:`k` value used to calculate the power spectrum (**default: max_k_calc=1e2**)
             - nk (*int*) : the number of :math:`k` values used to calculate the power spectrum (**default: nk=3000**)
             - h (*float*) : the dimensionless Hubble parameter (**default: h=0.6737**)
             - omega0_b (*float*) : the current baryon density, :math:`\Omega_{b,0}h^2` (**default: omega0_b=0.02233**)
             - omega0_cdm (*float*) : the current cold dark matter density, :math:`\Omega_{cdm,0}h^2` (**default: omega0_cdm=0.11977**)
             - n_s (*float*) : the tilt of the primordial power spectrum (**default: n_s=0.9652**)
             - sigma_8 (*float*) : the amplitude of matter fluctations within a sphere of radius :math:`r=8\ \mathrm{Mpc}\ \mathrm{h}^{-1}` at :math:`z=0` (**default: sigma_8=0.8101**)
             - verbose (*bool*) : whether to turn on the default CLASS logging (**default: verbose=False**)
             - gauge (*string*) : whether to use the synchronous or newtonian gauge (**default: gauge=sync**)
             - output (*string*) : whether to output the matter power spectrum or CMB power spectrum (**default: output='mPk'**)

.. py:function:: linear_growth_factor(z_val=0.0)

Function to calculate the linear growth factor :math:`D_1\left(z\right)` at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the linear growth factor is calculated (**default: z_val=0.0**)

.. py:function:: linear_power(input_k, z_val=0.0)

Function to calculate the linear power spectrum at a given redshift using Class

:Parameters: - input_k (*array*) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - z_val (*float*) : the redshift at which the linear power spectrum is calculated (**default: z_val=0.0**)

.. py:function:: zeldovich_power(n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False)

Function to compute the Zel'dovich approximation power spectrum at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the Zel'dovich power spectrum is calculated (**default: z_val=0.0**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=0`
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

.. important:: If input power spectrum values are given the log-spaced :math:`k` values used to compute it must also be given.

.. py:function:: lpt_one_loop_power(z_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False)

Function to calculate the LPT 1-loop power spectrum spectrum at a given redshift using `Wang et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.437..588W/abstract>`__

:Parameters: - z_val (*float*) : the redshift at which the Zel'dovich power spectrum is calculated (**default: z_val=0.0**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=0`
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

.. py:function:: gctm_power(n_val=32, zinit=99.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10))

Function to calculation the GCTM power spectrum at a given redshift using the Beyond Zel'dovich approximation if function :math:`A\left(z\right)` is not defined. See Lane et al. 2020 for more details

:Parameters: - z_val (*float*) : the redshift at which the GCTM power spectrum is calculated (**default: z_val=0.0**)
             - zinit (*float*) : the initial redshift :math:`z_i` that the time dependent functions are integrated from (**default: zinit=99.0**)
             - epsilon (*float*) : the value of the expansion parameter :math:`\epsilon` (**default: epsilon=1.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             -input_k_init (*array*, optional) : log-space array of :math:`k` values are which the input_P is calculated at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=z_i`
             - input_A (*array*, optional) : array of :math:`A\left(z\right)` values
             - input_z (*array*, optional) : array of redshift values corresponding to :math:`A` values
             - input_B (*array*, optional) : array of :math:`B\left(z\right)` values corresponding to :math:`A` values and redshift values
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

.. important:: If input_A or input_B or both are given then the redshift values used to compute them must also be given. For convergence use at least 1000 redshift values.

.. important:: If input_P is given you must also give input_k_init and if input_P is not evaluated at :math:`z_i=99` and you have not passed your own input_A array you must also specify the initial redshift at which input_P is calculated as z_init. 

Example I - Calculating the linear power spectrum
-------------------------------------------------

.. jupyter-execute::
    :hide-code:

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

.. jupyter-execute::

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

Example II - Calculating the Zel'dovich power spectrum
------------------------------------------------------

.. jupyter-execute::
    :hide-output:

    # Calculate the Zel'dovich power spectrum at z=0

    P_zel_0=GCTM(nk=1000).zeldovich_power(input_k=k_vals)

.. jupyter-execute::

  # Plot the results

  plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
  plt.loglog(k_vals, P_zel_0, color=colours[3], linestyle='--', linewidth=2.2, label=r"$\mathrm{Zel}^\prime\mathrm{dovich}$")
  plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
  plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
  plt.legend(loc="upper right", frameon=False, fontsize=14.)
  plt.xlim([1e-3, 1])
  plt.ylim([1e1, 1e5])
  plt.show()

We can also calculate the Zel'dovich power spectrum using a Gaussian damped initial power spectrum given by

.. math::

  \mathrm{P}_\mathrm{damped}\left(k\right)=\mathrm{e}^{-\left(\frac{k}{k_c}\right)^2}\mathrm{P}_\mathrm{lin}\left(k\right)

.. jupyter-execute::
    :hide-output:

    # Calculate the Zel'dovich power spectrum at z=0 with kc=5 h/Mpc

    P_zel_0_5=GCTM(nk=1000).zeldovich_power(input_k=k_vals, kc=5.0)

.. jupyter-execute::

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

Example III - Calculating the GCTM power spectrum
--------------------------------------------------

We can also calculate the Beyond Zel'dovich power spectrum if no :math:`A\left(z\right)` and :math:`B\left(z\right)` functions are specified. These functions are

.. math::

  A\left(z\right)=\frac{D_1\left(z\right)}{D_1\left(z_i\right)}\ \mathrm{and\ } B\left(z\right)=-\epsilon\frac{3}{2}H_0^2\Omega_m\int_z'^z\frac{dz''}{a''H\left(z''\right)}\int_{z_i}^z\frac{dz'}{H\left(z'\right)}\left(\frac{D_1\left(z\right)}{D_1\left(z_i\right)}\right)^2.

See Lane et al. (2020) for more details.

.. jupyter-execute::
    :hide-output:

    # Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc

    P_gctm_0_5=GCTM(nk=1000).gctm_power(input_k=k_vals, kc=5.0)

.. jupyter-execute::

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

You can also define your own :math:`A\left(z\right)` function. The :math:`B\left(z\right)` is calculated as

.. math::

  B\left(z\right)=-\epsilon\frac{3}{2}H_0^2\Omega_m\int_z'^z\frac{dz''}{a''H\left(z''\right)}\int_{z_i}^z\frac{dz'}{H\left(z'\right)}\left(A\left(z\right)\right)^2.

.. jupyter-execute::
    :hide-output:

    # Define redshift values

    z_vals=np.linspace(0.0, 200.0, 100)

    # Calculate A values

    A_vals=np.zeros_like(z_vals)

    for i in range(100):

      A_vals[i]=GCTM().linear_growth_factor(z_val=z_vals[i])/GCTM().linear_growth_factor(z_val=99.0)

    # Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc with input A

    P_gctm_input_A=GCTM(nk=1000).gctm_power(input_k=k_vals, kc=5.0, input_z=z_vals, input_A=A_vals)

.. jupyter-execute::

  # Plot the results

  plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
  plt.loglog(k_vals, P_gctm_input_A, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
  plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
  plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
  plt.legend(loc="lower left", frameon=False, fontsize=14.)
  plt.xlim([1e-3, 1])
  plt.ylim([1e1, 1e5])
  plt.show()
