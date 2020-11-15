import numpy as np
import sys
from scipy.interpolate import interp1d as interp
from scipy.integrate import trapz
import mcfit
from mcfit import SphericalBessel as sph
import warnings
import os
from GCTM.cosmology_tools.cosmology import Cosmo
from GCTM.cosmology_tools.time_dependence import TimeDep
from GCTM.power_spectrum_tools.integration import PowerIntegral
from GCTM.power_spectrum_tools.covariances import calc_covariances_func, calc_covariances_func_one_loop
from GCTM.misc_tools.timer import Timer
from GCTM.misc_tools.progress_update import progress_func_power

# Parameters: min_k is the minimum k value, max_k is the maximum k value,
# nk is the number of k values, h is the Hubble constant, omega0_b is the baryon density today,
# omega0_cdm is the CDM density today, n_s is the spectral index
# k_max is the maximum k value used when calculating the power spectrum


class RSD(object):

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

        # Define the k vector for integrals

        self.k_int = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        # Initialise Classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)


    # Function to calculate the Kaiser redshift space power spectrum P(k)=(1+f_1\mu_k**2)**2*P_L(k)

    def kaiser(self, k_vals, z_val=0.0, mu_k_val=0.0):

        kaiser = (1.0+self.cosmo.calc_independent_linear_growth(z_val)*mu_k_val**2)**2*self.cosmo.calc_linear_power(k_vals, z_val=z_val)

        return kaiser
