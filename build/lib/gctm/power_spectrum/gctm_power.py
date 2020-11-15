import numpy as np
import sys
from scipy.interpolate import interp1d as interp
import warnings
import os
from GCTM.cosmology_tools.cosmology import Cosmo
from GCTM.cosmology_tools.time_dependence import TimeDep
from GCTM.power_spectrum_tools.integration import PowerIntegral
from GCTM.power_spectrum_tools.covariances import calc_covariances_func
from GCTM.misc_tools.timer import Timer
from GCTM.misc_tools.progress_update import progress_func_power
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate

warnings.filterwarnings("ignore")


class PowerSpec:

    """

    Class to calculate the GCTM power spectrum

    Parameters:

    min_k = the minimum k value,
    max_k = the maximum k value,
    nk = the number of k values,
    h = H_0/100,
    omega0_b = $\Omega_bh^2$ the baryon density today,
    omega0_cdm = $\Omega_cdmh^2$ the CDM density today,
    n_s = the spectral index
    k_c = cutoff k value if using an initial Gaussian damped power spectrum
    input_k = input k values (if you do not specify then automatic k values are returned)
    input_P = input P values at z=0 you must also specify the k values used to calculate input_P
    z_val = the redshift value at which the GCTM power spectrum will be calculated
    epsilon = the expansion parameter (controls the size of the non-linear correction)
    z_init = the initial redshift value the time dependent factors are integrated from
    input_z = input z values
    input_A = input A values you must also specify the z values used to calculate input_A
    input_B = input B values you must also specify the z values used to calculate input_B

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

        self.k_int = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def calc_gctm_power(self, n_val=32., zinit=99.0, z_val=0.0, epsilon=1.0, save=False, kc=0.0, input_k=np.zeros(10), input_P=np.zeros(10), input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        # Function to calculate the GCTM power spectrum

        time = Timer()
        time.start()

        if input_k.all() != 0.0 and input_P.all() != 0.0:

            D_1 = self.cosmo.calc_linear_growth(zinit)

            if kc != 0.0:

                P_vals = np.exp(-(self.k_int/kc)**2)*(D_1)**2*input_P

            else:

                P_vals = (D_1)**2*input_P

            P_func = interp(input_k, P_vals)

        else:

            P = self.cosmo.calc_linear_power(self.k_int)

            D_1 = self.cosmo.calc_linear_growth(zinit)

            if kc != 0.0:

                P_vals = np.exp(-(self.k_int/kc)**2)*(D_1)**2*P

            else:

                P_vals = (D_1)**2*P

            P_func = interp(self.k_int, P_vals)

        print("Calculated the input power spectrum")

        if input_A.all() != 0.0 and input_B.all() == 0.0:

            A_func = interp(input_z, input_A)
            A_squared_val = A_func(z_val)**2

            B_squared_val = TimeDep().calc_B(input_A, z_val)

        if input_A.all() != 0.0 and input_B.all() != 0.0:

            A_func = interp(input_z, input_A)
            A_squared_val = A_func(z_val)**2

            B_func = interp(input_z, input_B)
            B_squared_val = B_func(z_val)**2

        if input_A.all() == 0.0 and input_B.all() == 0.0:

            D_1_val = self.cosmo.calc_linear_growth(z_val)
            D_1_init = self.cosmo.calc_linear_growth(zinit)

            A_squared_val = (D_1_val/D_1_init)**2

            w = 1.5*self.cosmo.Omega_m()*self.cosmo.H0()**2
            zeta_val = TimeDep().calc_zeta(z_val)
            B_squared_val = (epsilon*w*zeta_val)**2

        print('Calculated time dependent functions A(z) and B(z)')

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, sigma_0_0, D_vals, F_vals, G_vals, H_vals, I_vals, rho_vals = calc_covariances_func(self.k_int, P_func)

        print("Calculated the covariances")

        XY = X_vals + Y_vals

        W_prime = -(1.0/3.0)*eta_E*sigma_psi-(1.0/3.0)*(sigma_psi-0.5*X_vals)*(5.0*I_vals+8.0*H_vals+G_vals-(2.0/3.0)*sigma_0)

        Z_prime = (5.0*F_vals+D_vals)**2+(1.0/9.0)*rho_vals**2+(2.0/3.0)*rho_vals*(5.0*F_vals+D_vals)+(2.0/3.0)*(sigma_psi-0.5*X_vals)*(G_vals+7.0*H_vals)-(1.0/3.0)*Y_vals*(-2.5*I_vals+3.0*H_vals+0.5*G_vals+(1.0/3.0)*sigma_0)

        exponent_1 = A_squared_val*X_vals-2.0*B_squared_val*W_prime
        exponent_3 = A_squared_val*Y_vals-2.0*B_squared_val*Z_prime

        front = exponent_3
        exponent_k_squared = exponent_1+exponent_3
        zero_lag_1 = sigma_psi*(A_squared_val+(1.0/3.0)*B_squared_val*eta_E)

        P_calculated = np.zeros_like(self.k_int)

        for i in range(self.nk):

            progress_func_power(i, self.nk)

            P_calculated[i] = PowerIntegral().calc_power_rs(n_val, self.k_int[i], P_func, q_vals, A_squared_val, front, exponent_k_squared, zero_lag_1, sigma_psi, self.max_k)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(self.k_int, P_calculated)
            P_return = P_calculated_func(input_k)

            if save is True:

                for i in range(len(P_return)):

                    with open(os.path.join("results/P_gctm_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_return[i])+'\n')

            return P_return

        else:

            if save is True:

                for i in range(self.nk):

                    with open(os.path.join("results/P_gctm_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_calculated[i])+'\n')

            return self.k_int, P_calculated
