import numpy as np
import mcfit
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.special import gamma as Gamma_func
from scipy.special import hyp1f1
from scipy.special import hyperu
from scipy.special import factorial, factorial2


class ZeldovichPowerIntegral(mcfit.mcfit):

    """

    Class to calculate the Zel'dovich power spectrum using mcfit

    """

    def __init__(self, r, n):
        self.n = n
        UK = mcfit.kernels.Mellin_SphericalBesselJ(self.n)
        mcfit.mcfit.__init__(self, r, UK, q=1.5-n, lowring=True)

        # set pre and post factors
        self.prefac = r**(3-n)
        self.postfac = (2*np.pi)**1.5


class PowerIntegral:

    """

    Class to calculate the real and redshift space power spectra using mcfit

    Parameters:

    - n_val is the number of spherical-Bessel functions summed over
    - ki is the k value
    - P is the linear power spectrum function
    - q is the array of q values
    - A is the value of A(z)^2
    - front is the expression at the front of the integral
    - exponent_k_squared is the expression in the exponent multiplied by k^2
    - zero_lag_1 is the zero lag term
    - sigma_psi is the variance at zero displacement
    - kmax is the maximum k value to integrate the low-k approximation to

    """

    def calc_power_rs(self, n_val, ki, P, q, A, front, exponent_k_squared, zero_lag_1, sigma_psi, kmax):

        # Function to calculate the real space power spectrum

        Power = 0.0

        f = np.array([], dtype=np.float128)

        # Calculate the low-k approximation if ki<5e-5

        if ki < 5e-3:

            Power = calc_low_k_approx(ki, P, sigma_psi, A, kmax)

        else:

            for n in range(0, int(n_val)):

                I = ZeldovichPowerIntegral(q, n)

                if n > 0:

                    f = (ki*front)**n*np.exp(-0.5*ki**2*exponent_k_squared)

                else:

                    f = np.exp(-0.5*ki**2*exponent_k_squared) - np.exp(-ki**2*zero_lag_1)

                kk, this_Pzel = I(f, extrap=False)

                Power += spline(kk, this_Pzel)(ki)

        return Power

    def calc_power_rsd(self, ki, q, P, front, exponent_k_squared, zero_lag_1, kmax, alpha_0, alpha_1, C, xi, kappa, gamma, sigma, rsd_exponent_1, rsd_exponent_2, sigma_psi, A, f_val, mu_k):

        # Function to calculate the redshift space power spectrum

        Power = 0.0

        # Calculate the low k approximation if ki<5e-4

        if ki < 5e-5:

            Power = calc_low_k_approx_rsd(ki, P, sigma_psi, A, kmax, f_val, mu_k)

        else:
            for n in range(0, 32):

                K_n = 0.0
                I = ZeldovichPowerIntegral(q, n)

                if n > 0:

                    for l in range(0, 1):

                        K_n += np.exp(-0.5*ki**2*rsd_exponent_1-0.5*ki**2*rsd_exponent_2)*calc_F_l(l, alpha_0, alpha_1, -0.5*ki**2*C, xi, kappa, gamma, sigma)*(-1)**l*np.power(alpha_1, n)*np.power(kappa, n-l)*np.power(gamma**2/sigma, l)*hyperu(-l, n-l+1.0, 0.5*ki**2*C*kappa*alpha_1)


                    f = (ki*front)**n*np.exp(-0.5*ki**2*exponent_k_squared)*K_n

                else:

                    for l in range(0, 1):

                        K_n += np.exp(-0.5*ki**2*rsd_exponent_1-0.5*ki**2*rsd_exponent_2)*calc_F_l(l,alpha_0,alpha_1,-0.5*ki**2*C,xi,kappa,gamma,sigma)*(-1)**l*np.power(kappa,-l)*np.power(gamma**2/sigma,l)*hyperu(-l,-l+1.0,0.5*ki**2*C*kappa*alpha_1)


                    f = (np.exp(-0.5*ki**2*exponent_k_squared) - np.exp(-ki**2*zero_lag_1))*K_n

                kk, this_Pzel = I(f, extrap=False)

                Power += spline(kk, this_Pzel)(ki)

            return Power

def calc_low_k_approx(k, P, sigma_psi, A, kmax):

    # Low k approximation in real space

    Q3 = quad(lambda q: (P(q)/q)**2, 1e-5, kmax)[0]

    term_1 = (1.0-A*k**2*sigma_psi+0.5*A**2*k**4*sigma_psi**2)*A*P(k)
    term_2 = 0.5*(k**4)*Q3/(10.0*np.pi**2)

    return term_1+term_2


def calc_low_k_approx_rsd(k, P, sigma_psi, A, kmax, f_val, mu_k):

    # Low k approximation in redshift space

    factor = f_val*(f_val+2.0)*mu_k**2
    Q3 = quad(lambda q: (P(q)/q)**2, 1e-5, kmax)[0]

    term_1 = (1.0-A*k**2*sigma_psi-A*k**2*factor*sigma_psi+0.5*A**2*k**4*sigma_psi**2)*A*(1.0+f_val*mu_k**2)**2*P(k)+(1.0-A*k**2*sigma_psi-A*k**2*factor*sigma_psi)*((1.0+f_val*mu_k**2)**2*P(k)*k**2*(1.0+factor)*sigma_psi)
    term_2 = 0.5*(k**4)*Q3/(10.0*np.pi**2)

    return term_1+term_2


def calc_F_l(l, alpha_0, alpha_1, C, xi, kappa, gamma, sigma):

    # Calculate the function F_l for the redshift space power spectrum

    F_l = (np.power(np.pi, -0.5)*np.power(4.0, l)*Gamma_func(0.5))/(Gamma_func(1.0)*Gamma_func(1.0-l)*Gamma_func(2.0*l+1.0))*hyp1f1(l, l+0.5, C*alpha_1*(gamma**2/sigma))*hyp1f1(0.5, 1.0, C*alpha_0*(xi*gamma**2/(sigma*kappa)))

    for m in range(1, l):

        F_l += ((np.power(np.pi, -0.5)*np.power(-1.0, m)*np.power(4.0, l-m)*Gamma_func(m+0.5))/(Gamma_func(m+1.0)*Gamma_func(1.0+2.0*m-l)*Gamma_func(2.0*l-2.0*m+1.0)))*np.power(alpha_0/alpha_1, m)*np.power(xi/kappa, m)*hyp1f1(l-2.0*m, l-m+0.5, C*alpha_1*(gamma**2/sigma))*hyp1f1(m+0.5, m+1.0, C*alpha_0*(xi*gamma**2/(sigma*kappa)))

    return np.array(F_l, dtype=np.float128)
