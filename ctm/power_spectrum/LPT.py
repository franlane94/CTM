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
from .integration import PowerIntegral
from .covariances import calc_covariances_func
from ..misc.timer import Timer
from ..misc.progress_update import progress_func_power
from ..misc.progress_update import progress_func_power_save

warnings.filterwarnings("ignore")


class LPTPower:

    """

    Class to calculate the Zel'dovich power spectrum

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

        self.k = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        # Initialise Classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def calc_zeldovich_power(self, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        """
        Function to calculate the power spectrum for the Zel'dovich approximation using the method presented in 1209.0780 and 1410.1617
        """

        time = Timer()
        time.start()

        # Calculate the linear power spectrum for the specified redshift value

        if input_k.all() != 0.0 and input_P.all() != 0.0:

            D_1 = self.cosmo.calc_linear_growth(z_val)

            if kc != 0.0:

                P_vals = np.exp(-(input_k/kc)**2)*(D_1)**2*input_P

            else:

                P_vals = (D_1)**2*input_P

            P_func = interp(input_k, P_vals)

        else:

            P = self.cosmo.calc_linear_power(self.k, 0.0)

            D_1 = self.cosmo.calc_linear_growth(z_val)

            if kc != 0.0:

                P_vals = np.exp(-(self.k/kc)**2)*(D_1)**2*P

            else:

                P_vals = (D_1)**2*P

            P_func = interp(self.k, P_vals)

            k_int = self.k

            nk_calc = self.nk

            max_k_calc = self.max_k

        print("Calculated the input power spectrum")

        # Calculate the covariances

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, D_vals, F_vals, G_vals = calc_covariances_func(k_int, P_func)

        print("Calculated the covariances")

        XY = X_vals + Y_vals

        front = Y_vals
        exponent_k_squared = XY
        zero_lag_1 = sigma_psi

        A = 1.0

        # Begin calculating the Zel'dovich power spectrum

        P_calculated = np.zeros_like(k_int)

        for i in range(nk_calc):

            progress_func_power(i, nk_calc)

            P_calculated[i] = PowerIntegral().calc_power_rs(n_val, k_int[i], P_func, q_vals, A, front, exponent_k_squared, zero_lag_1, sigma_psi, max_k_calc)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(k_int, P_calculated)
            P_return = P_calculated_func(input_k)

            if save is True:

                for i in range(len(P_return)):

                    with open(os.path.join("P_zeldovich_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_return[i])+'\n')

            return P_return

        else:

            if save is True:

                for i in range(len(P_calculated)):

                    with open(os.path.join("P_zeldovich_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_calculated[i])+'\n')

            return k_int, P_calculated
