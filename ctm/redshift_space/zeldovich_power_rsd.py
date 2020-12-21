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
from ..power_spectrum.integration import PowerIntegral
from ..power_spectrum.covariances import calc_covariances_func, calc_covariances_func_one_loop
from ..misc.timer import Timer
from ..misc.progress_update import progress_func_power

warnings.filterwarnings("ignore")

# Parameters: min_k is the minimum k value, max_k is the maximum k value,
# nk is the number of k values, h is the Hubble constant, omega0_b is the baryon density today,
# omega0_cdm is the CDM density today, n_s is the spectral index
# k_max is the maximum k value used when calculating the power spectrum


class ZelPowerRSD(object):

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

    def calc_zeldovich_power_rsd(self, z_val=0.0,  mu_k_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the power spectrum for the Zel'dovich approximation in redshift space using the method presented in

        time = Timer()
        time.start()

        # Calculate the input linear power spectrum for the specified redshift

        if input_k.all() != 0.0 and input_P.all() != 0.0:

            D_1 = self.cosmo.calc_linear_growth(zinit)

            if kc != 0.0:

                P_vals = np.exp(-(input_k/kc)**2)*(D_1)**2*input_P

            else:

                P_vals = (D_1)**2*input_P

            P_func = interp(input_k, P_vals)

        else:

            P = self.cosmo.calc_linear_power(self.k_int, 0.0)

            D_1 = self.cosmo.calc_linear_growth(z_val)

            if kc != 0.0:

                P_vals = np.exp(-(self.k_int/kc)**2)*(D_1)**2*P

            else:

                P_vals = (D_1)**2*P

            P_func = interp(self.k_int, P_vals)

            k_calc = self.k_int

            nk_calc = self.nk

            max_k_calc = self.max_k

        print("Calculated the input power spectrum")

        # Calculate the covariances

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, D_vals, F_vals, G_vals = calc_covariances_func(k_calc, P_func)

        print("Calculated the covariances")

        # Define the arrays to be passed to the integral Class (see documentation for more details)

        XY = X_vals + Y_vals

        f_1_val = self.cosmo.calc_independent_linear_growth(z_val)

        C = Y_vals

        front = Y_vals
        exponent_k_squared = XY
        zero_lag_1 = sigma_psi

        rsd_exponent_1 = X_vals*f_1_val*mu_k_val**2*(2.0+f_1_val)

        rsd_exponent_2 = Y_vals*f_1_val*mu_k_val**2*(2.0+f_1_val*mu_k_val**2)

        # Define the time dependent factors A(z) and B(z) which is 1 and 0 respectively in the case of LPT

        A = 1.0

        # Define the additional RSD functions needed for GCTM which are all 1 in the case of LPT (see documentation for more details)

        xi = 1.0
        gamma = 1.0
        kappa = 1.0
        sigma = 1.0

        alpha_0, alpha_1, alpha_2, alpha_3 = calc_alpha_vals(self, mu_k_val, z_val)

        # Begin calculating the Zel'dovich power spectrum in redshift space

        P_calculated = np.zeros_like(k_calc)

        for i in range(nk_calc):

            progress_func_power(i, nk_calc)

            P_calculated[i] = PowerIntegral().calc_power_rsd(k_calc[i], q_vals, P_func, front, exponent_k_squared, zero_lag_1, self.max_k, alpha_0, alpha_1, C, xi, kappa, gamma, sigma, rsd_exponent_1, rsd_exponent_2, sigma_psi, A, f_1_val, mu_k_val)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(k_calc, P_calculated)
            P_return = P_calculated_func(input_k)

            if save is True:

                for i in range(len(P_return)):

                    with open(os.path.join("P_zeldovich_rsd_"+str(int(z_val))+"_"+str(int(mu_k_val))+".txt"), 'a') as file:

                        file.writelines(str(P_return[i])+'\n')

            return P_return

        else:

            if save is True:

                for i in range(len(P_calculated)):

                    with open(os.path.join("P_zeldovich_rsd_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_calculated[i])+'\n')

            return k_calc, P_calculated


# Function to calculate alpha_0, alpha_1, alpha_2 and alpha_3 as defined in the user documentation


def calc_alpha_vals(self, mu_k_val, z_val):

    alpha_0 = 1.0+self.cosmo.calc_independent_linear_growth(z_val)*mu_k_val**2*(2.0+self.cosmo.calc_independent_linear_growth(z_val))

    alpha_1 = (1.0+self.cosmo.calc_independent_linear_growth(z_val)*mu_k_val**2)**2

    alpha_2 = 2.0*self.cosmo.calc_independent_linear_growth(z_val)*mu_k_val*np.sqrt(1.0-mu_k_val**2)*(self.cosmo.calc_independent_linear_growth(z_val)*mu_k_val**2+1.0)

    alpha_3 = self.cosmo.calc_independent_linear_growth(z_val)**2*mu_k_val**2*(1.0-mu_k_val**2)

    return alpha_0, alpha_1, alpha_2, alpha_3
