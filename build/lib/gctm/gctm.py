import numpy as np
from cosmology_tools.cosmology import Cosmo
from power_spectrum_tools.LPT import LPTPower
from power_spectrum_tools.gctm_power import PowerSpec
from redshift_space_tools.zeldovich_power_rsd import ZelPowerRS
from redshift_space_tools.alternate_rsd import RSD
from redshift_space_tools.gctm_rsd import GCTMPowerRS


class GCTM:

    """

    Class to calculate the GCTM power spectrum

    Cosmological Parameters:
    min_k = the minimum k value,
    max_k = the maximum k value,
    nk = the number of k values,
    h = H_0/100,
    omega0_b = $\Omega_bh^2$ the baryon density today,
    omega0_cdm = $\Omega_cdmh^2$ the CDM density today,
    n_s = the spectral index

    Zel'dovich Parameters:

    input_k = input k values (if you do not specify then automatic k values are returned)
    input_P = input P values at z=0 you must also specify the k values used to calculate input_P
    z_val = the redshift value at which the Zel'dovich power spectrum will be calculated

    GCTM Parameters:

    k_c = cutoff k value if using an initial Gaussian damped power spectrum
    input_k = input k values (if you do not specify then automatic k values are returned)
    input_P = input P values at z=0 you must also specify the k values used to calculate input_P
    z_val = the redshift value at which the GCTM power spectrum will be calculated
    epsilon = the expansion parameter (controls the size of the non-linear correction)
    z_init = the initial redshift value the time dependent factors are integrated from
    input_z = input z values
    input_A = input A values you must also specify the z values used to calculate input_A
    input_B = input B values you must also specify the z values used to calculate input_B

    Redshift parameters:

    mu_k_val = the line-of-sight value i.e. \mu_k=\hat{k}_i\hat{z}_i where \hat{z}_i is the global line-of-sight

    """

    def __init__(self, min_k=1e-7, max_k=2e3, min_k_calc=1e-5, max_k_calc=100.0, nk=3000, h=0.6737, omega0_b=0.02242, omega0_cdm=0.11933, n_s=0.9665, sigma_8=0.8102, verbose=False, gauge='sync', output='mPk', **kwargs):

        self.min_k = min_k
        self.max_k = max_k
        self.min_k_zel = min_k_calc
        self.max_k_zel = max_k_calc
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

        # Function to calculate the linear growth factor using Classylss

        return Cosmo(h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, k_max=self.max_k_zel, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_linear_growth(z_val=z_val)

    def linear_power(self, input_k, z_val=0.0):

        # Function to calculate the linear power spectrum using Classylss

        return Cosmo(h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, k_max=self.max_k_zel, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_linear_power(input_k, z_val=z_val)

    def zeldovich_power(self, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the Zel'dovich power spectrum in real space

        # To save the output power spectrum set save=True

        return LPTPower(min_k=self.min_k, max_k=self.max_k, min_k_zel=self.min_k_zel, max_k_zel=self.max_k_zel, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_zeldovich_power(n_val=n_val, z_val=z_val, kc=kc, input_k=input_k, input_P=input_P, save=save)

    def lpt_one_loop_power(self, z_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the LPT 1-loop power spectrum in real space

        # To save the output power spectrum set save=True

        return LPTPower(min_k=self.min_k, max_k=self.max_k, min_k_zel=self.min_k_zel, max_k_zel=self.max_k_zel, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).one_loop(z_val=z_val, input_k=input_k, input_P=input_P, save=save)

    def gctm_power(self,  n_val=32, zinit=99.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):


        # Function to calculate the GCTM power spectrum in real space

        # If no A and B functions are specified the Beyond Zel'dovich approximation is used. See the Documentation for details

        return PowerSpec(min_k=self.min_k_zel, max_k=10.0, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_gctm_power(zinit=zinit, z_val=z_val, epsilon=epsilon, save=save, kc=kc, input_k=input_k, input_P=input_P, input_z=input_z, input_A=input_A, input_B=input_B)

    def zeldovich_power_rsd(self, z_val=0.0,  mu_k_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the Zel'dovich power spectrum in redshift space

        return ZelPowerRS(min_k=self.min_k_zel, max_k=10.0, nk=1500, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_zeldovich_power_rs(z_val=z_val,  mu_k_val=mu_k_val, input_k=input_k, input_P=input_P, save=save)

    def kaiser_power(self, input_k=np.zeros(10), z_val=0.0, mu_k_val=0.0):

        # Function to calculate the Kaiser power spectrum in redshift space

        return RSD(min_k=self.min_k_zel, max_k=10.0, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).kaiser(k_vals=input_k, z_val=z_val, mu_k_val=mu_k_val)

    def gctm_power_rsd(self, zinit=99.0, z_val=0.0, epsilon=1.0, kc=0.0, mu_k_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False, input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        # Function to calculate the GCTM power spectrum in redshift space

        # If no A and B functions are specified the Beyond Zel'dovich approximation is used. See Documentation page of the wiki for more information

        return GCTMPowerRS(min_k=self.min_k_zel, max_k=10.0, nk=1500, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_gctm_power_rs(zinit=zinit, z_val=z_val, epsilon=epsilon, kc=kc, mu_k_val=mu_k_val, input_k=input_k, input_P=input_P, save=save, input_z=input_z, input_A=input_A, input_B=input_B)
