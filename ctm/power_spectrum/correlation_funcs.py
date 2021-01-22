import numpy as np
import sys
from scipy.interpolate import interp1d as interp
from scipy.integrate import trapz
import mcfit
from mcfit import SphericalBessel as sph
import warnings
import os
from ..cosmology.cosmology import Cosmo
from ..cosmology.time_dependence import TimeDep
from .LPT import LPTPower
from .ctm_power import PowerSpec
from ..misc.timer import Timer
from ..misc.progress_update import progress_func_power
from ..misc.progress_update import progress_func_power_save

warnings.filterwarnings("ignore")

### Constant definitions

npi2 = np.power(np.pi,-2)
renorm = np.sqrt(0.5*np.pi)


class Corr:

    """

    Class to calculate the two-point correlation functions for the Zel'dovich approximation and second-order CTM

    Parameters:

    - min_k is the minimum k value
    - max_k is the maximum k value
    - nk is the number of k values
    - h is the Hubble constant
    - omega0_b is the baryon density today
    - omega0_cdm is the CDM density today
    - n_s is the spectral index
    - k_max is the maximum k value used when calculating the power spectrum

    - z_val is the redshift value at which the Zel'dovich power spectrum is calculated
    - input_k are user specified input k values (if you do not specify then automatic k values are returned)
    - input_P is an input power spectrum at z=0 if this is given then the k values used to calculate input_P must also be given as input_k
    - k_c is the cutoff k value if using an initial Gaussian damped power spectrum
    - n_val is the number of spherical Bessel functions summed over
    - epsilon is the expansion parameter which controls the size of the second-order CTM term
    - z_init is the initial redshift value the time dependent factors are integrated from
    - input_z is input z values
    - input_A is input A values you must also specify the z values used to calculate input_A
    - input_B is input B values you must also specify the z values used to calculate input_B

    - min_r is the minimum r value used in the calculation of the correlation function
    - max_r is the maximum r value used in the calculation of the correlation function
    - nr is the number of r values used in the calcultion of the correlation function

    """

    def __init__(self, min_k, max_k, nk, h, omega0_b, omega0_cdm, n_s, sigma_8, verbose, gauge, output, **kwargs):

        self.h = h
        self.omega0_b = omega0_b
        self.omega0_cdm = omega0_cdm
        self.n_s = n_s
        self.sigma_8 = sigma_8
        self.gauge = gauge
        self.output = output
        self.verbose = verbose
        self.nk = nk
        self.min_k = min_k
        self.max_k = max_k

        # Define the k vector for the integrals

        self.k_int = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        # Initialise classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def corr_func(self, k_values, P_values, min_r=1.0, max_r=1000.0, nr=10000):

        """
        Function to calculate the two-point correlation function given k_values and P_values
        """

        P_func = interp(k_values, P_values)

        new_k_vals = np.logspace(np.log10(1.0/max_r) , np.log10(1.0/min_r), nr)

        r_vals, I0 = dosph(0, new_k_vals, P_func, 0)

        corr_vals = 0.5*npi2*I0

        return r_vals, corr_vals

    def zeldovich_corr_func(self, min_r=1.0, max_r=1000.0, nr=10000, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        """
        Function to calculate the two-point correlation function for the Zel'dovich approximation
        """

        if input_k.all() == 0.0:

            k_values, P_values = LPTPower(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_zeldovich_power(n_val=n_val, z_val=z_val, kc=kc, input_k=input_k, input_P=input_P, save=save)

        else:

            k_values = input_k

            P_values = LPTPower(min_k=self.min_k, max_k=self.max_k, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_zeldovich_power(n_val=n_val, z_val=z_val, kc=kc, input_k=input_k, input_P=input_P, save=save)

        P_func = interp(k_values, P_values)

        new_k_vals = np.logspace(np.log10(1.0/max_r) , np.log10(1.0/min_r), nr)

        r_vals, I0 = dosph(0, new_k_vals, P_func, 0)

        corr_vals = 0.5*npi2*I0

        return r_vals, corr_vals

    def ctm_corr_func(self, min_r=1.0, max_r=1000.0, nr=10000, n_val=32, zinit=100.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_k_init=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        """
        Function to calculate the two-point correlation function for the CTM
        """

        if input_k.all() == 0.0:

            k_values, P_values = PowerSpec(min_k=self.min_k, max_k=10.0, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_ctm_power(n_val=n_val, zinit=zinit, z_val=z_val, epsilon=epsilon, save=save, kc=kc, input_k=input_k, input_P=input_P, input_k_init=input_k_init, input_z=input_z, input_A=input_A, input_B=input_B)

        else:

            k_values = input_k

            P_values = PowerSpec(min_k=self.min_k, max_k=10.0, nk=self.nk, h=self.h, omega0_b=self.omega0_b, omega0_cdm=self.omega0_cdm, n_s=self.n_s, sigma_8=self.sigma_8, verbose=self.verbose, gauge=self.gauge, output=self.output).calc_ctm_power(n_val=n_val, zinit=zinit, z_val=z_val, epsilon=epsilon, save=save, kc=kc, input_k=input_k, input_P=input_P, input_k_init=input_k_init, input_z=input_z, input_A=input_A, input_B=input_B)

        P_func = interp(k_values, P_values)

        new_k_vals = np.logspace(np.log10(1.0/max_r) , np.log10(1.0/min_r), nr)

        r_vals, I0 = dosph(0, new_k_vals, P_func, 0)

        corr_vals = 0.5*npi2*I0

        return r_vals, corr_vals

### Function to calculate spherical Bessel integrals returns q values and integral values
### Parameters: n = order of Bessel, x = variable of the function, f = function integrating over, a = power law

def dosph(n, x, f, a, tilt=1.5):

    func = renorm*np.power(x, a)*f(x)

    return sph(x, nu=n, q=tilt, lowring=True)(func, extrap=True)
