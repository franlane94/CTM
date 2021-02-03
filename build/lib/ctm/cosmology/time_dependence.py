import numpy as np
import sys
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.misc import derivative
from scipy.interpolate import interp1d as interp
from .cosmology import Cosmo
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
    - gauge = the perturbation gauge used
    - z_init = the initial redshift defined in the CTM

    - k = the k values to calculate the function at
    - z = the redshift value or values to calculate the time dependent functions at

    """

    def __init__(self, zinit=100.0, h=0.6737, omega0_b=0.02233, omega0_cdm=0.11933, k_max=10.0, n_s=0.9665, sigma_8=0.8102, verbose=False, gauge='sync', output='mPk', **kwargs):

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

        """
        Function to calculate \zeta (see documentation for more details)
        """

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

            def super_conformal(z_prime, z_val):

                return quad(lambda p: np.power(H(p), -1)*(D_1(p)/D_1(self.zinit))**2, self.zinit, z_prime)[0]

            for i in range(len(z)):

                zeta[i] = quad(lambda z_prime: super_conformal(z_prime, z[i])*np.power(scale_factor(z_prime)*H(z_prime),-1), self.zinit, z[i], epsabs=1e-4)[0]

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

                return quad(lambda p: np.power(H(p), -1)*(D_1(p)/D_1(self.zinit))**2, self.zinit, z_prime)[0]

            zeta = quad(lambda z_prime: super_conformal(z_prime, z)*np.power(scale_factor(z_prime)*H(z_prime),-1), self.zinit, z, epsabs=1e-4)[0]

        return zeta

    def calc_B(self, A_vals, z):

        """
        Function to calculate a general B(z) (see documentation for more details)
        """

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

            def super_conformal(z_prime, z_val):

                return quad(lambda p: np.power(H(p), -1)*A(p)**2, self.zinit, z_prime)[0]

            B = np.zeros_like(z)

            for i in range(int(len(z))):

                B[i] = quad(lambda z_prime: super_conformal(z_prime, z[i])*np.power(scale_factor(z_prime)*H(z_prime),-1), self.zinit, z[i], epsabs=1e-4)[0]

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

                return quad(lambda p: np.power(H(p), -1)*A(p)**2, self.zinit, z_prime)[0]

            B = quad(lambda z_prime: super_conformal(z_prime, z)*np.power(scale_factor(z_prime)*H(z_prime),-1), self.zinit, z, epsabs=1e-4)[0]

        return B
