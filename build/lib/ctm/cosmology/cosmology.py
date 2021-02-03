import numpy as np
import sys
import classylss
import classylss.binding as CLASS


class Cosmo:

    """
    Class to calculate cosmological functions

    - k_max = the maximum k value used when calculating the power spectra
    - h = H_0/100
    - omega0_b = $\Omega_bh^2$ the baryon density today
    - omega0_cdm = $\Omega_cdmh^2$ the CDM density today
    - n_s = the spectral index
    - gauge = the perturbation gauge used

    - z_val = the redshift value or values to calculate the function at
    - k = the k values to calculate the function at

    """

    def __init__(self, h, omega0_b, omega0_cdm, k_max, n_s, sigma_8, verbose, gauge, output, **kwargs):

        # Read in the cosmological parameters as defined in ctm.py

        self.h = h
        self.omega0_b = omega0_b
        self.omega0_cdm = omega0_cdm
        self.P_k_max = k_max+1.0
        self.n_s = n_s
        self.sigma_8 = sigma_8
        self.gauge = gauge
        self.output = output

        # Initialise Classylss

        self.engine = CLASS.ClassEngine({'h': self.h, 'omega_cdm': self.omega0_cdm, 'omega_b': self.omega0_b, 'n_s': self.n_s, 'gauge': self.gauge,
        'output': self.output, 'sigma8': self.sigma_8, 'P_k_max_1/Mpc': self.P_k_max})

    def scale_factor(self, z_val):

        """
        Function to calculate the scale factor a=(1+z)^{-1}
        """

        return np.power(1.0+z_val, -1)

    def calc_hubble(self, z_val):

        """
        Function to calculate the Hubble function H(z)
        """

        H0 = self.h*100.0

        bg = CLASS.Background(self.engine)

        H = H0*bg.efunc(z_val)

        return H

    def calc_linear_growth(self, z_val):

        """
        Function to calculate the linear growth factor D_1
        """

        bg = CLASS.Background(self.engine)

        D_1 = bg.scale_independent_growth_factor(z_val)

        return D_1

    def calc_independent_linear_growth(self, z_val):

        """
        Function to calculated f=dlnD_1/dlna
        """

        bg = CLASS.Background(self.engine)

        f = bg.scale_independent_growth_rate(z_val)

        return f

    def calc_linear_power(self, k,  z_val=0.0):

        """
        Function to calculate the linear power spectrum
        """

        bg = CLASS.Background(self.engine)

        sp = CLASS.Spectra(self.engine)

        D_1 = bg.scale_independent_growth_factor(z_val)

        pk_lin = D_1**2*sp.get_pklin(k, 0.0)

        return pk_lin

    def Omega_m(self):

        """
        Function to calculate \Omega_m
        """

        return (self.omega0_cdm+self.omega0_b)/self.h**2

    def H0(self):

        """
        Function to calculate H0
        """

        return self.h*100.0
