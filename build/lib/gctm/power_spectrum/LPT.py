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
from GCTM.misc_tools.progress_update import progress_func_power_save

warnings.filterwarnings("ignore")


class LPTPower:

    """

    Class the calculate the Zel'dovich and LPT 1-loop power spectra

    Parameters:

    - min_k is the minimum k value
    - max_k is the maximum k value
    - nk is the number of k values
    - h is the Hubble constant
    - omega0_b is the baryon density today
    - omega0_cdm is the CDM density today
    - n_s is the spectral index
    - k_max is the maximum k value used when calculating the power spectrum


    """

    def __init__(self, min_k, max_k, min_k_zel, max_k_zel, nk, h, omega0_b, omega0_cdm, n_s, sigma_8, verbose, gauge, output, **kwargs):

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
        self.min_k_zel = min_k_zel
        self.max_k_zel = max_k_zel

        # Define the k vector for loop integrals

        self.k = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        self.k_zel = np.logspace(np.log10(self.min_k_zel), np.log10(self.max_k_zel), self.nk)

        # Initialise Classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def calc_zeldovich_power(self, n_val=32, z_val=0.0, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the power spectrum for the Zel'dovich approximation using the method presented in 1209.0780 and 1410.1617

        time = Timer()
        time.start()

        # Calculate the linear power spectrum for the specified redshift value

        if input_k.all() != 0.0 and input_P.all() != 0.0:

            D_1 = self.cosmo.calc_linear_growth(zinit)

            if kc != 0.0:

                P_vals = np.exp(-(input_k/kc)**2)*(D_1)**2*input_P

            else:

                P_vals = (D_1)**2*input_P

            P_func = interp(input_k, P_vals)

        else:

            P = self.cosmo.calc_linear_power(self.k_zel, 0.0)

            D_1 = self.cosmo.calc_linear_growth(z_val)

            if kc != 0.0:

                P_vals = np.exp(-(self.k_zel/kc)**2)*(D_1)**2*P

            else:

                P_vals = (D_1)**2*P

            P_func = interp(self.k_zel, P_vals)

            k_calc = self.k_zel

            nk_calc = self.nk

            max_k_calc = self.max_k_zel

        print("Calculated the input power spectrum")

        # Calculate the covariances

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, sigma_0_0, D_vals, F_vals, G_vals, H_vals, I_vals, rho_vals = calc_covariances_func(k_calc, P_func)

        print("Calculated the covariances")

        XY = X_vals + Y_vals

        front = Y_vals
        exponent_k_squared = XY
        zero_lag_1 = sigma_psi

        A = 1.0

        # Begin calculating the Zel'dovich power spectrum

        P_calculated = np.zeros_like(k_calc)

        for i in range(nk_calc):

            progress_func_power(i, nk_calc)

            P_calculated[i] = PowerIntegral().calc_power_rs(n_val, k_calc[i], P_func, q_vals, A, front, exponent_k_squared, zero_lag_1, sigma_psi, max_k_calc)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(k_calc, P_calculated)
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

            return k_calc, P_calculated

    def one_loop(self, z_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False):

        # Function to calculate the power spectrum for LPT+1-loop using the method presented in 1209.0780 and 1410.1617

        # Load weight parameters for the Q and R integrals
        full_path = os.path.realpath(__file__)
        glxval, glwval = np.loadtxt(os.path.dirname(full_path)+"/gl_128.txt", unpack=True)
        glxval = glxval.reshape(1, -1)
        glwval = glwval.reshape(1, -1)

        # Calculate the linear power spectrum for the loop integrals

        P_lin = self.cosmo.calc_linear_power(self.k)
        P_int = interp(self.k, P_lin)

        time = Timer()
        time.start()

        if input_k.all() != 0.0 and input_P.all() != 0.0:

            P_vals = input_P

            P_func = interp(input_k, P_vals)

            k_calc = input_k

            nk_calc = len(input_k)

            max_k_calc = max(input_k)

        else:

            P = self.cosmo.calc_linear_power(self.k_zel, 0.0)

            D_1 = self.cosmo.calc_linear_growth(z_val)

            P_vals = (D_1)**2*P

            P_func = interp(self.k_zel, P_vals)

            k_calc = self.k_zel

            nk_calc = self.nk

            max_k_calc = self.max_k_zel

        print("Calculated the input power spectrum")

        # Calculate the Q and R functions as defined in 1209.0780

        k_loop = np.logspace(-4, np.log10(50), self.nk)

        Q = np.zeros_like(k_loop)
        R = np.zeros_like(k_loop)

        for i in range(len(k_calc)):

            Q[i] = Q_external(k_loop[i], P_int, glxval, glwval)
            R[i] = R_external(k_loop[i], P_int, glxval, glwval)

        Q_func = interp(k_loop, Q)
        R_func = interp(k_loop, R)

        sigma_psi, sigma_22, sigma_13, q_vals, X_loop, Y_loop = calc_covariances_func_one_loop(k_loop, P_func, Q_func, R_func)

        print("Calculated the covariances")

        XY = X_loop + Y_loop

        # Define the arrays to be passed to the integral Class (see documentation for more details)

        front = Y_loop
        exponent_k_squared = XY
        zero_lag_1 = sigma_psi + sigma_22 + 2.0*sigma_13

        # Define the time dependent factor A(z) which is 1 in the case of LPT

        A = 1.0

        # Begin calculating the LPT 1-loop power spectrum

        P_calculated = np.zeros_like(k_calc)

        for i in range(nk_calc):

            progress_func_power(i, nk_calc)

            P_calculated[i] = PowerIntegral().calc_power_rs(k_calc[i], P_func, q_vals, A, front, exponent_k_squared, zero_lag_1, sigma_psi, max_k_calc)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(k_calc, P_calculated)
            P_return = P_calculated_func(input_k)

            if save is True:

                for i in range(len(P_return)):

                    with open(os.path.join("P_lpt_oneloop_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_return[i])+'\n')

            return P_return

        else:

            if save is True:

                for i in range(len(P_calculated)):

                    with open(os.path.join("P_lpt_oneloop_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_calculated[i])+'\n')

            return self.k, P_calculated


def Q1(r, x):

    # Function to calculate the Q1 function defined in equation (B17) in Arxiv:1209.0780 for LPT 1-loop

    y = 1.0+r**2-2.0*r*x
    Q_1 = r**2*(1.0-x**2)**2/y**2
    return Q_1


def Q_internal(k, r, P, glxval, glwval):

    # Function to calculate the internal integral in equation (B17) in Arxiv:1209.0780 for LPT 1-loop

    def func(r, x):

        return P(k*np.sqrt(1.0+r**2-2.0*r*x))*Q1(r, x)

    return (glwval*func(r, glxval)).sum(axis=-1)


def Q_external(k, P, glxval, glwval):

    # Function to calculate the external integral in equation (B17) in Arxiv:1209.0780 for LPT 1-loop

    fac = k**3/(2.0*np.pi)**2
    kint = np.logspace(-5, 2, 1000)
    r = (kint/k)
    tol = 1e-5
    absr1 = abs(r-1)
    mask = absr1 < tol

    y = np.zeros_like(r)
    y[~mask] = Q_internal(k, r[~mask].reshape(-1, 1), P, glxval, glwval)
    y *= P(kint)

    return fac*trapz(y, r)


def R1(r, x):

    # Function to calculate the R1 function defined in equation (B16) in Arxiv:1209.0780 for LPT 1-loop

    y = 1.0+r**2-2.0*r*x

    return (r**2*(1.0-x**2)**2)/y


def R_internal(r, P, glxval, glwval):

    # Function to calculate the internal integral defined in equation (B16) in Arxiv:1209.0780 for LPT 1-loop

    def func(r, x):

        return R1(r, x)

    return (glwval*func(r, glxval)).sum(axis=-1)


def R_external(k, P, glxval, glwval):

    # Function to calculate the external integral defined in equation (B16) in Arxiv:1209.0780 for LPT 1-loop

    fac = k**3/(2.0*np.pi)**2*P(k)
    kint = np.logspace(-4, np.log10(5e2), 1000)
    r = (kint/k)
    absr1 = abs(r-1)
    tol = 1e-5
    mask = absr1 < tol

    y = np.zeros_like(r)
    y[~mask] = R_internal(r[~mask].reshape(-1, 1), P, glxval, glwval)
    y *= P(kint)

    return fac*trapz(y, r)
