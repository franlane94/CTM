import numpy as np
import sys
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.misc import derivative
from scipy.interpolate import interp1d as interp
from GCTM.cosmology_tools.cosmology import Cosmo
import warnings
warnings.filterwarnings("ignore")


class TimeDep:

    """
    Class to calculate the time dependent functions

    Parameters:

    - h is the Hubble constant
    - omega0_b is the baryon fraction today
    - omega0_cdm is the CDM fraction today
    - k_max is the maximum k value used when calculating the power spectra
    - n_s is the spectral index

    """

    def __init__(self, zinit=99.0, h=0.6737, omega0_b=0.02233, omega0_cdm=0.11977, k_max=10.0, n_s=0.9652, sigma_8=0.8101, verbose=False, gauge='sync', output='mPk', **kwargs):

        # Define the initial z value

        self.zinit = zinit

        # Read in the cosmological parameters as specified in gctm.py

        self.h = h
        self.omega0_b = omega0_b
        self.omega0_cdm = omega0_cdm
        self.k_max = k_max+1.0
        self.n_s = n_s
        self.sigma_8 = sigma_8
        self.gauge = gauge
        self.output = output
        self.verbose = verbose

        # Initialise classylss

        self.cosmo = Cosmo(self.h, self.omega0_b, self.omega0_cdm, self.k_max, self.n_s, self.sigma_8, self.verbose, self.gauge, self.output)

    def calc_zeta(self, z):

        # Function to calculate \zeta (see documentation for more details)

        if type(z) == np.ndarray:

            # Calculate the Hubble parameter for the redshifts

            hubble_values = self.cosmo.calc_hubble(z)
            H = interp(z, hubble_values)

            # Calculate the linear growth factor D_1 values for the redshifts

            linear_growth_values = self.cosmo.calc_linear_growth(z)
            D_1 = interp(z, linear_growth_values)

            # Calculate the scale factor values for the redshifts

            scale_factor_values = self.cosmo.scale_factor(z)
            scale_factor = interp(z, scale_factor_values)

            # Calculate \zeta

            zeta = np.zeros_like(z)

            for i in range(int(len(z))):

                zeta[i] = dblquad(lambda s,p: np.power(H(s),-1)*np.power(scale_factor(p)*H(p),-1)*(D_1(p)/D_1(self.zinit))**2, self.zinit, z[i], lambda p: p, lambda p: z[i], epsabs=1e-4)[0]

        else:

            # Define the redshifts

            z_vals = np.linspace(0.0, 200.0, 1000)

            # Calculate the Hubble parameter for the redshifts

            hubble_values = self.cosmo.calc_hubble(z_vals)
            H = interp(z_vals, hubble_values)

            # Calculate the linear growth factor D_1 values for the redshifts

            linear_growth_values = self.cosmo.calc_linear_growth(z_vals)
            D_1 = interp(z_vals, linear_growth_values)

            # Calculate the scale factor values for the redshifts

            scale_factor = interp(z_vals, self.cosmo.scale_factor(z_vals))

            # Calculate \zeta

            def super_conformal(z_prime, z_val):

                return quad(lambda p: np.power(scale_factor(p)*H(p),-1), z_prime, z_val)[0]

            zeta = quad(lambda z_prime: super_conformal(z_prime, z)*np.power(H(z_prime),-1)*(D_1(z_prime)/D_1(self.zinit))**2, self.zinit, z, epsabs=1e-4)[0]

        return zeta

    def calc_B(self, A_vals, z):

        # Function to calculate a general B(z) (see documentation for more details)

        if type(z) == np.ndarray:

            # Calculate the Hubble parameter for the redshifts

            hubble_values = self.cosmo.calc_hubble(z)
            H = interp(z, hubble_values)

            # Calculate the linear growth function D_1 for the redshifts

            linear_growth_values = self.cosmo.calc_linear_growth(z)
            A = interp(z, A_vals)

            # Calculate the scale factor values for the redshifts

            scale_factor_values = self.cosmo.scale_factor(z)
            scale_factor = interp(z, scale_factor_values)

            # Calculate B

            def super_conformal(z_prime, z):

                return quad(lambda p: np.power(scale_factor(p)*H(p),-1), z_prime, z)[0]

            B = np.zeros_like(z)

            for i in range(int(len(z))):

                B[i] = quad(lambda z_prime: super_conformal(z_prime, z[i])*np.power(H(z_prime),-1)*A(z_prime)**2, self.zinit, z[i], epsabs=1e-4)[0]

        else:

            # Define the redshifts

            z_vals = np.linspace(0.0, 200.0, len(A_vals))

            # Calculate the Hubble parameter for the redshifts

            hubble_values = self.cosmo.calc_hubble(z_vals)
            H = interp(z_vals, hubble_values)

            # Calculate the linear growth function D_1 for the redshifts

            linear_growth_values = self.cosmo.calc_linear_growth(z_vals)

            # Interpolate the array of A(z) values

            A = interp(z_vals, A_vals)

            # Calculate the scale factor for the redshifts

            scale_factor = interp(z_vals, self.cosmo.scale_factor(z_vals))

            # Calculate B(z)

            def super_conformal(z_prime, z_val):

                return quad(lambda p: np.power(scale_factor(p)*H(p),-1), z_prime, z_val)[0]

            B = quad(lambda z_prime: super_conformal(z_prime, z)*np.power(H(z_prime),-1)*A(z_prime)**2, self.zinit, z, epsabs=1e-4)[0]

        return B

"""

    def calc_alpha(self, z):

        # Function to calculate \alpha (see documentation for more details)

        self.nz = int(len(z))

        hubble_values = self.cosmo.calc_hubble(z)
        H = interp(z, hubble_values)

        linear_growth_values = self.cosmo.calc_linear_growth(z)
        D_1 = interp(z, linear_growth_values)

        D1deriv = np.diff(linear_growth_values)/np.diff(z)
        D1deriv_func = interp(z[0:self.nz-1], D1deriv)

        int3 = np.zeros_like(z)

        for i in range(self.nz):

            int3[i] = quad(lambda t: np.power(H(t), -1), self.zinit, z[i], limit=1000)[0]

        alpha = np.zeros_like(z, dtype=np.float128)

        for i in range(self.nz):

            alpha[i] = 1.0+H(self.zinit)*(D1deriv_func(self.zinit)/D_1(self.zinit))*int3[i]

        return alpha

    def calc_beta(self, z):

        # Function to calculate \beta (see documentation for more details)

        self.nz = int(len(z))

        hubble_values = cosmo.calc_hubble(z)
        H = interp(z, hubble_values)

        linear_growth_values = cosmo.calc_linear_growth(z)
        D_1 = interp(z, linear_growth_values)

        beta = np.zeros_like(z, dtype=np.float128)

        for i in range(self.nz):

            beta[i] = dblquad(lambda s, p: np.power(H(s), -1)*(D_1(p)/D_1(self.zinit))*np.power(H(p)*cosmo.scale_factor(p), -1), self.zinit, z[i], lambda p: p, lambda p:z, epsabs=1e-4)[0]

        return beta

    def calc_gamma(self,z):

        # Function to calculate \gamma (see documentation for more details)

        self.nz = int(len(z))

        if self.nz < 1000:

            return print('For convergence please use 1000 or more redshift values')

        hubble_values = cosmo.calc_hubble(z)
        H = interp(z, hubble_values)

        linear_growth_values = cosmo.calc_linear_growth(z)
        D_1 = interp(z, linear_growth_values)

        gamma = np.zeros_like(z,dtype=np.float128)

        for i in range(self.nz):

            gamma[i] = dblquad(lambda s,p: np.power(H(s),-1)*(D_1(p)/D_1(self.zinit))*np.power(H(p)*cosmo.scale_factor(p),-1)*calc_alpha(p), self.zinit, z, lambda p: p, lambda p: z, epsabs=1e-4)[0]

        return gamma
"""
