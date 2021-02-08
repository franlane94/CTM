import numpy as np
from .cosmology.cosmology import Cosmo
from .power_spectrum.LPT import LPTPower
from .power_spectrum.ctm_power import PowerSpec
from .power_spectrum.correlation_funcs import Corr

class CTM:

    """

    Class to calculate the CTM power spectrum

    Cosmological Parameters:

    - min_k = the minimum k value,
    - max_k = the maximum k value,
    - nk = the number of k values,
    - h = H_0/100,
    - omega0_b = $\Omega_bh^2$ the baryon density today,
    - omega0_cdm = $\Omega_cdmh^2$ the CDM density today,
    - n_s = the spectral index

    Zel'dovich Parameters:

    - input_k = input k values (if you do not specify then automatic k values are returned)
    - input_P = input P values at z=0 you must also specify the k values used to calculate input_P
    - z_val = the redshift value at which the Zel'dovich power spectrum will be calculated
    - k_c = cutoff k value if using an initial Gaussian damped power spectrum
    - n_val = the number of spherical Bessel functions summed over (n_val=32.)

    CTM Parameters:

    - k_c = cutoff k value if using an initial Gaussian damped power spectrum
    - input_k = input k values (if you do not specify then automatic k values are returned)
    - input_P = input P values at zinit you must also specify the k values used to calculate input_P as input_k_init. Note if your initial z is not z_init=100 then you must also specify z_init unless you have passed input_A
    - input_k_init = input k values at which the input_P power spectrum is calculated at
    - z_val = the redshift value at which the CTM power spectrum will be calculated
    - epsilon = the expansion parameter which controls the size of the second-order CTM terms (epsilon=1)
    - z_init = the initial redshift value the time dependent factors are integrated from (z_init=100)
    - input_z = input z values
    - input_A = input A values you must also specify the z values used to calculate input_A
    - input_B = input B values you must also specify the z values used to calculate input_B

    Correlation function parameters:

    - min_r = the minimum r value used to calculate the correlation function
    - max_r = the maximum r value used to calculate the correlation function
    - nr = the number of r values used to calculate the correlation function

    """

    def __init__(self, min_k=1e-5, max_k=100.0, nk=3000, h=0.6737, omega0_b=0.02233, omega0_cdm=0.11933, n_s=0.9665, sigma_8=0.8102, verbose=False, gauge='sync', output='mPk', **kwargs):

        self.min_k = min_k
        self.max_k = max_k
        self.nk = nk
        self.h = h
        self.omega0_b = omega0_b
        self.omega0_cdm = omega0_cdm
        self.n_s = n_s
        self.sigma_8 = sigma_8
        self.verbose = verbose
        self.gauge = gauge
        self.output = output

    def linear_growth_factor(self, z_val=0.0):

        """
        Function to calculate the linear growth factor using Classylss
        """

        return Cosmo(h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, k_max=self.max_k, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_linear_growth(z_val=z_val)

    def linear_power(self, input_k, z_val=0.0):

        """
        Function to calculate the linear power spectrum using Classylss
        """

        return Cosmo(h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, k_max=self.max_k, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_linear_power(input_k, z_val=z_val)

    def zeldovich_power(self, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        """
        Function to calculate the Zel'dovich power spectrum in real space.
        """

        return LPTPower(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_zeldovich_power(n_val=n_val, z_val=z_val, kc=kc, input_k=input_k, input_P=input_P, save=save)

    def ctm_power(self,  n_val=32, zinit=100.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_k_init=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        """
        Function to calculate the CTM power spectrum in real space.

        If no A and B functions are specified the Beyond Zel'dovich approximation is used. See the Documentation for details.
        """
        return PowerSpec(min_k=self.min_k, max_k=10.0, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_ctm_power(n_val=n_val, zinit=zinit, z_val=z_val, epsilon=epsilon, save=save, kc=kc, input_k=input_k, input_P=input_P, input_k_init=input_k_init, input_z=input_z, input_A=input_A, input_B=input_B)

    def corr_func(self, k_values, P_values, min_r=1.0, max_r=1000.0, nr=10000):

        """
        Function to calculate the two-point correlation function given k values and power spectrum values
        """

        return Corr(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).corr_func(min_r=min_r, max_r=max_r, nr=nr, k_values=k_values, P_values=P_values)

    def corr_func_zel(self, min_r=1.0, max_r=1000.0, nr=10000, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        """
        Function to calculate the Zel'dovich two-point correlation function in real space.
        """

        return Corr(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).zeldovich_corr_func(min_r=min_r, max_r=max_r, nr=nr, z_val=z_val, kc=kc, input_k=input_k, input_P=input_P, save=save)

    def corr_func_ctm(self, min_r=1.0, max_r=1000.0, nr=10000, n_val=32, zinit=100.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_k_init=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        """
        Function to calculate the CTM two-point correlation function in real space

        If no A and B functions are specified the Beyond Zel'dovich approximation is used. See the Documentation for details.
        """

        return Corr(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).ctm_corr_func(min_r=min_r, max_r=max_r, nr=nr, n_val=n_val, zinit=zinit, z_val=z_val, epsilon=epsilon, save=save, kc=kc, input_k=input_k, input_P=input_P, input_k_init=input_k_init, input_z=input_z, input_A=input_A, input_B=input_B)
