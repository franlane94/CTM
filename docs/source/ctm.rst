The CTM in real space
=====================
---------------------------------------------------------------------------------------------------------------


A brief introduction to the CTM and the Beyond Zel'dovich approximation
-----------------------------------------------------------------------

The CTM trajectory to second-order is given by

.. math::

  \mathbf{x}\left(\mathbf{q},t\right)=\mathbf{q}+\Psi^{\left(0\right)}_i\left(\mathbf{q},t_i\right)\left[A\left(t\right)\delta_{ij}+B_\epsilon\left(t\right)\bar{E}_{ij}\left(\mathbf{q},t_i\right)\right]

where the tidal tensor is :math:`\bar{E}_{ij}\left(\mathbf{q},t\right)=\nabla_i\nabla_j\nabla^{-2}\delta^{\left(0\right)}\left(\mathbf{q},t\right)` with :math:`\delta^{\left(0\right)}` being the linear overdensity field. This tidal term describes the effects of gravitational scattering and deflection. The time-dependent function :math:`B_\epsilon\left(t\right)` is,

.. math::

  B_\epsilon\left(t\right)=-\epsilon_\mathrm{BZ}\omega_0^2\int_{t_i}^{t}\frac{dt'}{a'^2}\int_{t_i}^{t'}\frac{dt''}{a''}A\left(t''\right)^2,

which after using :math:`dt=-\frac{a}{H}dz` can be written as,

.. math::

    B_\epsilon\left(t\right)=-\epsilon_\mathrm{BZ}\omega_0^2\int_{z_i}^{z}\frac{dz'}{a'H\left(z'\right)}\int_{z_i}^{z'}\frac{dz''}{H\left(z''\right)}A\left(z''\right)^2.

The Beyond Zel'dovich approximation is the second-order CTM trajectory given in with the following time-dependent functions

.. math::

  A\left(z\right)=\frac{D_1\left(z\right)}{D_1\left(z_i\right)}\ \mathrm{and\ }
    B_\epsilon\left(z\right)&=-\epsilon_\mathrm{BZ}\omega_0^2\int_{z_i}^{z}\frac{dz'}{a'H\left(z'\right)}\int_{z_i}^{z'}\frac{dz''}{H\left(z''\right)}\left(\frac{D_1\left(z''\right)}{D_1\left(z_i\right)}\right)^2.


Using the CTM module in real space
----------------------------------

Documentation introducing the CTM class.

.. py:class:: CTM(min_k=1e-5, max_k=100.0, nk=3000, h=0.6737, omega0_b=0.02233, omega0_cdm=0.11933, n_s=0.9665, sigma_8=0.8102, verbose=False, gauge='sync', output='mPk', **kwargs)

Class to calculate the linear, Zel'dovich and CTM power spectra in real space

:Parameters: - min_k (*int*) : the minimum :math:`k` value used to calculate the power spectrum (**default:       min_k=1e-5*)
             - max_k (*int*) : the maximum :math:`k` value used to calculate the power spectrum (**default: max_k=1e2**)
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

:Returns: - linear_growth_val (*float*) : the linear growth factor value :math:`D_1` at the given redshift

.. py:function:: linear_power(input_k, z_val=0.0)

Function to calculate the linear power spectrum at a given redshift using ``Classylss``

:Parameters: - input_k (*array*) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - z_val (*float*) : the redshift at which the linear power spectrum is calculated (**default: z_val=0.0**)

:Returns: - linear_power_spec (*array*) : array of :math:`P\left(k\right)` values

.. py:function:: zeldovich_power(n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False)

Function to compute the Zel'dovich approximation power spectrum at a given redshift

:Parameters: - z_val (*float*) : the redshift at which the Zel'dovich power spectrum is calculated (**default: z_val=0.0**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=0`
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

:Returns: - zel_power (*array*) : array of :math:`P\left(k\right)` values
          - k_values (*array*) : an array of log-spaced `k` values if input_k is not given

.. important:: If input power spectrum values are given the log-spaced :math:`k` values used to compute it must also be given.

.. py:function:: ctm_power(n_val=32, zinit=100.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_k_init=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10))

Function to calculation the second-order CTM power spectrum at a given redshift using the Beyond Zel'dovich approximation if function :math:`A\left(z\right)` is not defined. See Lane et al. 2021 for more details.

:Parameters: - z_val (*float*) : the redshift at which the CTM power spectrum is calculated (**default: z_val=0.0**)
             - zinit (*float*) : the initial redshift :math:`z_i` that the time dependent functions are integrated from (**default: zinit=100.0**)
             - epsilon (*float*) : the value of the expansion parameter :math:`\epsilon_\mathrm{BZ}` (**default: epsilon=1.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_k_init (*array*, optional) : log-space array of :math:`k` values are which the input_P is calculated at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=z_i`
             - input_A (*array*, optional) : array of :math:`A\left(z\right)` values
             - input_z (*array*, optional) : array of redshift values corresponding to :math:`A` values
             - input_B (*array*, optional) : array of :math:`B\left(z\right)` values corresponding to :math:`A` values and redshift values
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

:Returns: - zel_power (*array*) : array of :math:`P\left(k\right)` values
          - k_values (*array*) : an array of log-spaced `k` values if input_k is not given

.. important:: If input_A or input_B or both are given then the redshift values used to compute them must also be given. For convergence use at least 1000 redshift values.

.. important:: If input_P is given you must also give input_k_init and if input_P is not evaluated at :math:`z_i=100` and you have not passed your own input_A array you must also specify the initial redshift at which input_P is calculated as z_init.

.. py:function:: corr_func(k_values, P_values, min_r=1.0, max_r=1000.0, nr=10000)

Function to calculate the two-point correlation function given k values and a power spectrum

:Parameters: - k_values (*array*) : log-spaced array of :math:`k` values to compute the two-point correlation function at
            - P_values (*array*) : array of :math:`P\left(k\right)` values to compute the two-point correlation function with
            - min_r (*float*) : the minimum :math:`r` value returned (**default: min_r=1.0**)
            - max_r (*float*) : the maximum :math:`r` value returned (**default: min_r=1000.0**)
            - nr (*int*) : the number of :math:`r` values returned (**default: nr=10000**)

:Returns: - r_values (*array*) : array of :math:`r` values
          - corr_values (*array*) : array of :math:`\xi\left(r\right)` values

.. py:function:: corr_func_zel(min_r=1.0, max_r=1000.0, nr=10000, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False)

Function to calculate the two-point correlation function for the Zel'dovich approximation

:Parameters: - min_r (*float*) : the minimum :math:`r` value returned (**default: min_r=1.0**)
             - max_r (*float*) : the maximum :math:`r` value returned (**default: min_r=1000.0**)
             - nr (*int*) : the number of :math:`r` values returned (**default: nr=10000**)
             - z_val (*float*) : the redshift at which the Zel'dovich correlation function is calculated (**default: z_val=0.0**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=0`
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

:Returns: - r_values (*array*) : array of :math:`r` values
          - corr_values (*array*) : array of :math:`\xi\left(r\right)` values

.. important:: If input power spectrum values are given the log-spaced :math:`k` values used to compute it must also be given.

.. py:function:: corr_func_ctm(self, min_r=1.0, max_r=1000.0, nr=10000, n_val=32, zinit=100.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_k_init=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10))

Function to calculate the two-point correlation function for the CTM

:Parameters: - min_r (*float*) : the minimum :math:`r` value returned (**default: min_r=1.0**)
             - max_r (*float*) : the maximum :math:`r` value returned (**default: min_r=1000.0**)
             - nr (*int*) : the number of :math:`r` values returned (**default: nr=10000**)
             - z_val (*float*) : the redshift at which the CTM correlation function is calculated (**default: z_val=0.0**)
             - zinit (*float*) : the initial redshift :math:`z_i` that the time dependent functions are integrated from (**default: zinit=100.0**)
             - epsilon (*float*) : the value of the expansion parameter :math:`\epsilon_\mathrm{BZ}` (**default: epsilon=1.0**)
             - n_val (*int*, optional) : the number of spherical Bessel functions summed over (**default: n_val=32**)
             - kc (*float*, optional) : the cutoff value to be used if an initial Gaussian damped power spectrum is used (**default: kc=0.0**)
             - input_k (*array*, optional) : log-spaced array of :math:`k` values to compute the linear power spectrum at
             - input_k_init (*array*, optional) : log-space array of :math:`k` values are which the input_P is calculated at
             - input_P (*array*, optional) : array of linear power spectrum values at :math:`z=z_i`
             - input_A (*array*, optional) : array of :math:`A\left(z\right)` values
             - input_z (*array*, optional) : array of redshift values corresponding to :math:`A` values
             - input_B (*array*, optional) : array of :math:`B\left(z\right)` values corresponding to :math:`A` values and redshift values
             - save (*bool*, optional) : whether to save the calculated power spectrum to the current directory

:Returns: - r_values (*array*) : array of :math:`r` values
          - corr_values (*array*) : array of :math:`\xi\left(r\right)` values

.. important:: If input_A or input_B or both are given then the redshift values used to compute them must also be given. For convergence use at least 1000 redshift values.

.. important:: If input_P is given you must also give input_k_init and if input_P is not evaluated at :math:`z_i=100` and you have not passed your own input_A array you must also specify the initial redshift at which input_P is calculated as z_init.

Examples
--------

Example I - Calculating the linear power spectrum
*************************************************

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

Example II - Calculating the Zel'dovich power spectrum
******************************************************


.. jupyter-execute::
    :hide-output:

    # Calculate the Zel'dovich power spectrum at z=0

    P_zel_0=CTM(nk=1000).zeldovich_power(input_k=k_vals)

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

    P_zel_0_5=CTM(nk=1000).zeldovich_power(input_k=k_vals, kc=5.0)

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

Example III - Calculating the CTM power spectrum
************************************************

We can also calculate the Beyond Zel'dovich power spectrum if no :math:`A\left(z\right)` and :math:`B\left(z\right)` functions are specified. These functions are

.. math::

  A\left(z\right)=\frac{D_1\left(z\right)}{D_1\left(z_i\right)}\ \mathrm{and\ } B\left(z\right)=-\epsilon\frac{3}{2}H_0^2\Omega_m\int_z'^z\frac{dz''}{a''H\left(z''\right)}\int_{z_i}^z\frac{dz'}{H\left(z'\right)}\left(\frac{D_1\left(z\right)}{D_1\left(z_i\right)}\right)^2.

See Lane et al. (2021) for more details.

.. jupyter-execute::
    :hide-output:

    # Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc

    P_ctm_0_5=CTM(nk=1000).ctm_power(input_k=k_vals, kc=5.0)

.. jupyter-execute::

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

      A_vals[i]=CTM().linear_growth_factor(z_val=z_vals[i])/CTM().linear_growth_factor(z_val=99.0)

    # Calculate the Beyond Zel'dovich power spectrum at z=0 with kc=5 h/Mpc with input A

    P_ctm_input_A=CTM(nk=1000).ctm_power(input_k=k_vals, kc=5.0, input_z=z_vals, input_A=A_vals)

.. jupyter-execute::

  # Plot the results

  plt.loglog(k_vals, P_lin_0, color="black", linestyle='-', linewidth=2.2, alpha=0.8, label=r"$\mathrm{Linear}$")
  plt.loglog(k_vals, P_ctm_input_A, color=colours[2], linestyle='-.', linewidth=2.2, label=r"$\mathrm{Damped\ Beyond\ Zel}^\prime\mathrm{dovich}$")
  plt.xlabel(r"$k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]$", fontsize=14.)
  plt.ylabel(r"$\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]$", fontsize=14.)
  plt.legend(loc="lower left", frameon=False, fontsize=14.)
  plt.xlim([1e-3, 1])
  plt.ylim([1e1, 1e5])
  plt.show()

Example IV - Computing two-point correlation functions
******************************************************

.. jupyter-execute::
  :hide-output:

  # Compute the linear correlation function
  r_lin, corr_lin=CTM().corr_func(k_vals, P_lin_0)

  # Compute the Zel'dovich and CTM correlation functions

  r_zel, corr_zel=CTM(nk=1000).corr_func_zel()
  r_ctm, corr_ctm=CTM(nk=1000).corr_func_ctm()

.. jupyter-execute::

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
