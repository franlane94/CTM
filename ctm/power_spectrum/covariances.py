import numpy as np
import mcfit
from mcfit import SphericalBessel as sph
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.special import spherical_jn
from scipy.interpolate import interp1d as interp
from scipy.misc import derivative
import sys


# Constant definitions

npi2 = np.power(np.pi, -2)
renorm = np.sqrt(0.5*np.pi)


def calc_covariances_func(k, P):

    """
    Function to calculate the covariances X, Y, D, F, G, \sigma_psi,
    \sigma_0 and \eta_E for an input power spectrum
    """

    sigma_psi = calc_sigma_psi(k, P)
    q_vals, X_vals, Y_vals = calc_X_and_Y(k, P, sigma_psi)
    eta_E = calc_eta(k, P)
    _, D_vals, F_vals, sigma_0 = calc_D_and_F(k, P)
    _, G_vals = calc_G(k, P)

    return sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, D_vals, F_vals, G_vals

# Function definitions


def dosph(n, x, f, a, tilt=1.5):

    """
    Function to calculate spherical Bessel integrals returns q values and integral values
    Parameters: n = order of Bessel, x = variable of the function, f = function integrating over, a = power law
    """

    func = renorm*np.power(x, a)*f(x)

    return sph(x, nu=n, q=tilt, lowring=True)(func, extrap=True)


def calc_sigma_psi(k, P):

    """
    Function to calculate the variance of the displacement field \sigma^2{\psi} = \frac{1}{6\pi^2}\int{dkP_L\left(k)}
    """

    q, I0 = dosph(0, k, P, -2)
    sigma_vals = npi2*I0/6.0

    sigma_func = interp(q, sigma_vals)
    sigma_0 = sigma_func(min(q))

    return sigma_0


def calc_X_and_Y(k, P, sigma):

    """
    Function to calculate X and Y where X = \frac{1}{2\pi^2}\int{dk[\frac{2}{3}-\frac{j_1(kq)}{kq}]P_L} and Y = \frac{1}{2\pi^2}\int{dk[6\frac{j_1(kq)}{kq}-j_0(kq)]P_L} and sigma is \sigma^2_{\psi}
    """

    q, I0 = dosph(0, k, P, -2)
    _, I2 = dosph(2, k, P, -2)

    X = 2.0*sigma - npi2*(1.0/3.0)*(I0+I2)
    Y = npi2*I2

    return q, X, Y


def calc_eta(k, P):

    """
    Function to calculate eta^2_E = \frac{1}{6\pi^2}\int{dk}k^2P_L
    """

    q, I0 = dosph(0, k, P, 0)
    eta_vals = npi2*I0/6.0

    eta_func = interp(q, eta_vals)
    eta_0 = eta_func(min(q))

    return eta_0

def calc_D_and_F(k, P):

    """
    Function to calculate D=\frac{1}{2\pi^2}\int{dk kj_3(kq)P_L} and F=-\frac{1}{2\pi^2}\int{dk kj_2(kq)/kqP_L} and sigma^2_0=\frac{1}{\left(2\pi\right)^3}\int{dk}e^{-ik\cdotq}P_L
    """

    q, I0 = dosph(0, k, P, 0)
    _, I2 = dosph(2, k, P, 0)

    F = -0.5*npi2*I2
    D = npi2*(I0+I2)/(6.0)
    sigma_0 = 0.5*npi2*I0

    return q, D, F, sigma_0

def calc_G(k, P):

    """
    Function to calculate G=-\frac{1}{2\pi^2}\int{dk kj_1(kq)kP_L}
    """

    q, I1 = dosph(1, k, P, -1)

    G = -0.5*npi2*I1

    return q, G
