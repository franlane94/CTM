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
from ..power_spectrum.covariances import calc_covariances_func
from ..misc.timer import Timer
from ..misc.progress_update import progress_func_power
warnings.filterwarnings("ignore")


class CTMPowerRSD:

    """

    Class to calculate the CTM power spectrum

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

        # Define the k vector for integrals

        self.k_int = np.logspace(np.log10(self.min_k), np.log10(self.max_k), self.nk)

        # Initialise Classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def calc_ctm_power_rsd(self, zinit=100.0, z_val=0.0, epsilon=1.0, kc=0.0, mu_k_val=0.0, input_k=np.zeros(10), input_P=np.zeros(10), save=False, input_z=np.zeros(10), input_A=np.zeros(10), input_B=np.zeros(10)):

        time = Timer()
        time.start()

        if input_P.all() != 0.0:

            if kc != 0.0:

                P_vals = np.exp(-(self.k_int/kc)**2)*input_P

            else:

                P_vals = input_P

            P_func = interp(input_k_init, P_vals)

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

            B_squared_val = TimeDep(zinit, self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output).calc_B(input_A, z_val)

            z_vals = input_z

        if input_A.all() != 0.0 and input_B.all() != 0.0:

            A_func = interp(input_z, input_A)
            A_squared_val = A_func(z_val)**2

            B_func = interp(input_z, input_B)
            B_squared_val = B_func(z_val)**2

            z_vals = input_z

        if input_A.all() == 0.0 and input_B.all() == 0.0:

            D_1_val = self.cosmo.calc_linear_growth(z_val)
            D_1_init = self.cosmo.calc_linear_growth(zinit)

            A_squared_val = (D_1_val/D_1_init)**2

            w = 1.5*self.cosmo.Omega_m()*self.cosmo.H0()**2
            zeta_val = TimeDep(zinit, self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output).calc_zeta(z_val)

            B_squared_val = (epsilon*w*zeta_val)**2

            z_vals = np.linspace(0.0, 200.0, 100)

            B_array = TimeDep(zinit, self.h, self.omega0_b, self.omega0_cdm, self.max_k+1.0, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output).calc_zeta(z_vals)

        print('Calculated time dependent functions A(z) and B(z)')

        # Calculate the covariances

        sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, D_vals, F_vals, G_vals = calc_covariances_func(self.k_int, P_func)

        print("Calculated the covariances")

        # Define the arrays to be passed to the integral Class (see documentation for more details)

        XY = X_vals + Y_vals

        W_prime = -(1.0/3.0)*eta_E*sigma_psi+(sigma_psi-0.5*X_vals)*((2.0/3.0)*D_vals-(1.0/9.0)*sigma_0)

        Z_prime = (1.0/9.0)*G_vals**2+(2.0/3.0)*F_vals*(sigma_psi-0.5*X_vals)+Y_vals*((1.0/18.0)*sigma_0-(1.0/3.0)*D_vals-(1.0/3.0)*F_vals)

        f_1_val = self.cosmo.calc_independent_linear_growth(z_val)/D_1_init

        f_2_val = -epsilon*w*calc_derivative(z_vals, len(B_array), B_array, z_val)

        alpha_0, alpha_1, alpha_2, alpha_3 = calc_alpha_vals(f_1_val, mu_k_val)
        alpha_0_1, alpha_1_1, alpha_2_1, alpha_3_1 = calc_alpha_vals(f_2_val, mu_k_val)

        C = A_squared_val*Y_vals

        if alpha_2==0 and alpha_3==0:

            print(mu_k_val, alpha_2, alpha_3)
            xi = 1.0-2.0*(B_squared_val/A_squared_val)*(W_prime/X_vals)*(alpha_0_1/alpha_0)

            kappa = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_1_1/alpha_1)

            gamma = 1.0

            sigma = 1.0

        else:

            xi = 1.0-2.0*(B_squared_val/A_squared_val)*(W_prime/X_vals)*(alpha_0_1/alpha_0)

            kappa = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_1_1/alpha_1)

            gamma = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_2_1/alpha_2)

            sigma = 1.0-2.0*(B_squared_val/A_squared_val)*(Z_prime/Y_vals)*(alpha_3_1/alpha_3)

        exponent_1 = A_squared_val*X_vals-2.0*B_squared_val*W_prime
        exponent_3 = A_squared_val*Y_vals-2.0*B_squared_val*Z_prime

        front = exponent_3
        exponent_k_squared = exponent_1+exponent_3
        zero_lag_1 = sigma_psi*(A_squared_val+(1.0/3.0)*B_squared_val*eta_E)

        rsd_exponent_1 = A_squared_val*X_vals*f_1_val*mu_k_val**2*(2.0+f_1_val)-2.0*B_squared_val*W_prime*f_2_val*mu_k_val**2*(2.0+f_2_val)

        rsd_exponent_2 = A_squared_val*Y_vals*f_1_val*mu_k_val**2*(2.0+f_1_val*mu_k_val**2)-2.0*B_squared_val*Z_prime*f_2_val*mu_k_val**2*(2.0+f_2_val*mu_k_val**2)

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

def calc_derivative(z_vals, nz, array, z_val):

	a_vals = np.power(1.0+z_vals, -1)
	deriv = np.diff(np.log(array))/np.diff(np.log(a_vals))
	deriv_func = interp(a_vals[0:nz-1], deriv)

	a_vals_new = a_vals[0:nz-1]
	z_vals_new = z_vals[0:nz-1]
	calc_deriv = interp(z_vals_new,deriv_func(a_vals_new))

	return calc_deriv(z_val)

### Function to calculate alpha_0 and alpha_1

def calc_alpha_vals(f_1_val, mu_k_val):

	alpha_0 = 1.0+f_1_val*mu_k_val**2*(2.0+f_1_val)

	alpha_1 = (1.0+f_1_val*mu_k_val**2)**2

	alpha_2 = 2.0*f_1_val*mu_k_val*np.sqrt(1.0-mu_k_val**2)*(f_1_val*mu_k_val**2+1.0)

	alpha_3 = f_1_val**2*mu_k_val**2*(1.0-mu_k_val**2)

	return alpha_0, alpha_1, alpha_2, alpha_3

### Function to calculate alpha^1_0 and alpha^1_1

def calc_alpha_vals_1(f_2_func, mu_k_val):

	alpha_0 = 1.0+f_2_val*mu_k_val**2*(2.0+f_2_val)

	alpha_1 = (1.0+fun_2_val*mu_k_val**2)**2

	alpha_2 = 2.0*fun_2_val*mu_k_val*np.sqrt(1.0-mu_k_val**2)*(fun_2_val*mu_k_val**2+1.0)

	alpha_3 = fun_2_val**2*mu_k_val**2*(1.0-mu_k_val**2)

	return alpha_0, alpha_1, alpha_2, alpha_3
