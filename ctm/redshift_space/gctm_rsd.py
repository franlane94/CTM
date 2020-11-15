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


warnings.filterwarnings("ignore")

# Parameters: min_k is the minimum k value, max_k is the maximum k value,
# nk is the number of k values, h is the Hubble constant, omega0_b is the baryon density today,
# omega0_cdm is the CDM density today, n_s is the spectral index
# k_max is the maximum k value used when calculating the power spectrum


class GCTMPowerRS(object):

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

    def calc_gctm_power_rs(self, zinit=99.0, z_val=0.0, epsilon=1.0, kc=0.0, mu_k_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False, input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        time = Timer()
        time.start()

        # Calculate the input linear power spectrum for the specified redshift

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

        # Calculate the time dependent functions in A(z) and B(z)

        if input_A.all() != 0.0 and input_B.all() == 0.0:

            A_func = interp(input_z, input_A)
            A_squared_val = A_func(z_val)**2

            w = 1.5*self.cosmo.Omega_m()*self.cosmo.H0()**2
            B_squared_val = (epsilon*w*TimeDep().calc_B(input_A, z_val))**2

            z_vals = input_z

            A_array = input_A
            B_array = epsilon*w*TimeDep().calc_B(input_A, z_vals)

        if input_A.all() != 0.0 and input_B.all() != 0.0:

            A_func = interp(input_z, input_A)
            A_squared_val = A_func(z_val)**2

            w = 1.5*self.cosmo.Omega_m()*self.cosmo.H0()**2
            B_func = interp(input_z, input_B)
            B_squared_val = (epsilon*w*B_func(z_val))**2

            z_vals = input_z

            A_array = input_A
            B_array = epsilon*w*input_B

        if input_A.all() == 0.0 and input_B.all() == 0.0:

            D_1_val = self.cosmo.calc_linear_growth(z_val)
            D_1_init = self.cosmo.calc_linear_growth(zinit)

            A_squared_val = (D_1_val/D_1_init)**2

            w = 1.5*self.cosmo.Omega_m()*self.cosmo.H0()**2
            zeta_val = TimeDep().calc_zeta(z_val)
            B_squared_val = (epsilon*w*zeta_val)**2

            z_vals = np.linspace(0.0, 200.0, 1000)

            A_array = self.cosmo.calc_linear_growth(z_vals)/self.cosmo.calc_linear_growth(zinit)
            B_array = epsilon*w*TimeDep().calc_zeta(z_vals)

        print('Calculated time dependent functions A(z) and B(z)')

        # Calculate the covariances

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, sigma_0_0, D_vals, F_vals, G_vals, H_vals, I_vals, rho_vals = calc_covariances_func(self.k_int, P_func)

        print("Calculated the covariances")

        # Define the arrays to be passed to the integral Class (see documentation for more details)

        XY = X_vals + Y_vals

        W_prime = -(1.0/3.0)*eta_E*sigma_psi-(1.0/3.0)*(sigma_psi-0.5*X_vals)*(5.0*I_vals+8.0*H_vals+G_vals-(2.0/3.0)*sigma_0)

        Z_prime = (5.0*F_vals+D_vals)**2+(1.0/9.0)*rho_vals**2+(2.0/3.0)*rho_vals*(5.0*F_vals+D_vals)+(2.0/3.0)*(sigma_psi-0.5*X_vals)*(G_vals+7.0*H_vals)-(1.0/3.0)*Y_vals*(-2.5*I_vals+3.0*H_vals+0.5*G_vals+(1.0/3.0)*sigma_0)

        f_1 = calc_derivative(z_vals, len(A_array), A_array)
        f_2 = calc_derivative(z_vals, len(B_array), B_array)

        f_1_val = f_1(z_val)
        f_2_val = f_2(z_val)

        if mu_k_val==0.0:

            mu_k_val=1e-11

        if mu_k_val==1.0:

            mu_k_val=1.0-1e-11

        alpha_0, alpha_1, alpha_2, alpha_3 = calc_alpha_vals(f_1 , mu_k_val, z_val)
        alpha_0_1, alpha_1_1, alpha_2_1, alpha_3_1 = calc_alpha_vals(f_2, mu_k_val, z_val)

        C = A_squared_val*Y_vals

        xi = 1.0-2.0*(B_squared_val/A_squared_val)*(W_prime/X_vals)*(alpha_0_1/alpha_0)

        kappa = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_1_1/alpha_1)

        gamma = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_2_1/alpha_2)

        sigma = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_3_1/alpha_3)

        if np.isnan(gamma[0])==True:

            gamma = np.ones_like(X_vals)

        if np.isnan(sigma[0])==True:

            sigma = np.ones_like(X_vals)

        exponent_1 = A_squared_val*X_vals-2.0*B_squared_val*W_prime
        exponent_3 = A_squared_val*Y_vals-2.0*B_squared_val*Z_prime

        front = exponent_3
        exponent_k_squared = exponent_1+exponent_3
        zero_lag_1 = sigma_psi*(A_squared_val+(1.0/3.0)*B_squared_val*eta_E)

        rsd_exponent_1 = A_squared_val*X_vals*(f_1_val*mu_k_val**2*(2.0+f_1_val)-2.0*(B_squared_val/A_squared_val)*(W_prime/X_vals)*f_2_val*mu_k_val**2*(2.0+f_2_val))

        rsd_exponent_2 = A_squared_val*Y_vals*(f_1_val*mu_k_val**2*(2.0+f_1_val*mu_k_val**2)-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*f_2_val*mu_k_val**2*(2.0+f_2_val*mu_k_val**2))

        # Begin calculating the GCTM power spectrum in redshift space

        P_calculated = np.zeros_like(self.k_int)

        for i in range(self.nk):

            progress_func_power(i, self.nk)

            P_calculated[i] = PowerIntegral().calc_power_rsd(self.k_int[i], q_vals, P_func, front, exponent_k_squared, zero_lag_1, self.max_k, alpha_0, alpha_1, C, xi, kappa, gamma, sigma, rsd_exponent_1, rsd_exponent_2, sigma_psi, A_squared_val, f_1_val, mu_k_val)

        time.stop()

        if input_k.all() != 0.0:

            P_calculated_func = interp(self.k_int, P_calculated)
            P_return = P_calculated_func(input_k)

            if save is True:

                for i in range(len(P_return)):

                    with open(os.path.join("P_gctm_rsd_"+str(int(z_val))+"_"+str(int(mu_k_val))+".txt"), 'a') as file:

                        file.writelines(str(P_return[i])+'\n')

            return P_return

        else:

            if save is True:

                for i in range(len(P_calculated)):

                    with open(os.path.join("P_gctm_rsd_"+str(int(z_val))+".txt"), 'a') as file:

                        file.writelines(str(P_calculated[i])+'\n')

            return self.k, P_calculated

### Function to calculate derivatives dlnA/dlna

def calc_derivative(z_vals, nz, array):

	a_vals = np.power(1.0+z_vals, -1)
	deriv = np.diff(np.log(array))/np.diff(np.log(a_vals))
	deriv_func = interp(a_vals[0:nz-1], deriv)

	a_vals_new = a_vals[0:nz-1]
	z_vals_new = z_vals[0:nz-1]
	calc_deriv = interp(z_vals_new,deriv_func(a_vals_new))

	return calc_deriv

### Function to calculate alpha_0 and alpha_1

def calc_alpha_vals(f_1_func, mu_k_val, z_val):

	alpha_0 = 1.0+f_1_func(z_val)*mu_k_val**2*(2.0+f_1_func(z_val))

	alpha_1 = (1.0+f_1_func(z_val)*mu_k_val**2)**2

	alpha_2 = 2.0*f_1_func(z_val)*mu_k_val*np.sqrt(1.0-mu_k_val**2)*(f_1_func(z_val)*mu_k_val**2+1.0)

	alpha_3 = f_1_func(z_val)**2*mu_k_val**2*(1.0-mu_k_val**2)

	return alpha_0, alpha_1, alpha_2, alpha_3

### Function to calculate alpha^1_0 and alpha^1_1

def calc_alpha_vals_1(f_2_func, mu_k_val, z_val):

	alpha_0 = 1.0+f_2_func(z_val)*mu_k_val**2*(2.0+f_2_func(z_val))

	alpha_1 = (1.0+f_2_func(z_val)*mu_k_val**2)**2

	alpha_2 = 2.0*f_2_func(z_val)*mu_k_val*np.sqrt(1.0-mu_k_val**2)*(f_2_func(z_val)*mu_k_val**2+1.0)

	alpha_3 = f_2_func(z_val)**2*mu_k_val**2*(1.0-mu_k_val**2)

	return alpha_0, alpha_1, alpha_2, alpha_3
