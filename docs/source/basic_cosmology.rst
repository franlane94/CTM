Basic cosmology operations
=============================

Documentation introducing the Cosmo class

.. py:class:: Cosmo(h, omega0_b, omega0_cdm, k_max, n_s, sigma_8, verbose, gauge, output, **kwargs)

Class to calculate basic cosmology functions

:Parameters: - h (*float*) : the dimensionless Hubble parameter
             - omega0_b (*float*) : the current baryon density, :math:`\Omega_{b,0}h^2`
             - omega0_cdm (*float*) : the current cold dark matter density, :math:`\Omega_{cdm,0}h^2`
             - k_max (*float*) : the maximum :math:`k` value used when calculating the linear power spectrum
             - n_s (*float*) : the tilt of the primordial power spectrum
             - sigma_8 (*float*) : the amplitude of matter fluctuations within a sphere of radius :math:`r=8\ \mathrm{Mpc}\ \mathrm{h}^{-1}` at :math:`z=0`
             - verbose (*bool*) : whether to turn on the default ``class`` logging
             - gauge (*string*) : whether to use the synchronous or newtonian gauge
             - output (*string*) : whether to output the matter power spectrum or CMB power spectrum

.. py:function:: scale_factor(z_val)

Function to calculate the scale factor defined as :math:`a=\left(1+z\right)^{-1}`

:Parameters: - z_val (*float*) : the redshift at which the scale factor is computed

:Returns: - a_val (*float*) : the scale factor value at z_val

.. py:function:: calc_hubble(z_val)

Function to calculate the Hubble parameter, :math:`H\left(z\right)`, at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the Hubble parameter is computed

:Returns: - hubble_val (*float*) : the Hubble parameter at z_val

.. py:function:: Omega_m()

Function to calculate :math:`\Omega_m=\Omega_{cdm}+\Omega_b`

:Returns: - omega_m (*float*) : the :math:`\Omega_m` value

.. py:function:: H0

Function to calculate :math:`H_0=h\times100`

:Returns: - H0 (*float*) : the :math:`H_0` value

.. py:function:: calc_linear_growth(z_val)

Function to calculate the linear growth factor, :math:`D_1\left(z\right)`, at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the linear growth factor is computed

:Returns: - linear_growth_val (*float*) : the :math:`D_1` value at z_val

.. py:function:: calc_independent_linear_growth(z_val)

Function to calculate the independent linear growth factor, :math:`f=\frac{d\ln{D_1}}{d\ln{a}}`, at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the linear growth factor is computed

:Returns: - independent_growth_val (*float*) : the :math:`f` value at z_val

.. py:function:: calc_linear_power(k, z_val)

Function to calculate the linear power spectrum at a given redshift using ``classylss``

:Parameters: - z_val (*float*) : the redshift at which the linear growth factor is computed
             - k (*array*, *float*) : log-spaced array of :math:`k` values

:Returns: - linear_power (*array*) : the linear power spectrum evaluated at the given :math:`k` values and redshift value

Examples
--------

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

    linestyles = ["-", "--", "-.", ":"]

We can plot the linear growth factor for a range of redshifts

.. jupyter-execute::

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

We can also plot the independent growth factor for a range of redshifts

.. jupyter-execute::

  # Calculate the independent growth factor values

  f_vals=Cosmo(h=0.7, omega0_b=0.02233, omega0_cdm=0.112, n_s=0.96, sigma_8=0.8, k_max=10.0, verbose=False, gauge='sync', output='mPk').calc_independent_linear_growth(z_vals)

  # Plot the results

  plt.plot(z_vals, f_vals, color="black", linestyle='-', linewidth=2.2, alpha=0.8)
  plt.xlabel(r"$z$", fontsize=14.)
  plt.ylabel(r"$f$", fontsize=14.)
  plt.xlim([-2, 202])
  plt.ylim([0.5, 1.2])
  plt.show()
