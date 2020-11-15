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

    # Function to calculate the covariances X, Y, U, V, R, S, T
    # for an input power spectrum

    sigma_psi = calc_sigma_psi(k, P)
    q_vals, X_vals, Y_vals = calc_X_and_Y(k, P, sigma_psi)
    eta_E = calc_eta(k, P)
    sigma_0, sigma_0_0 = calc_sigma_0(k, P)
    _, D_vals, F_vals = calc_D_and_F(k, P)
    _, G_vals, H_vals, I_vals = calc_G_H_and_I(k, P)
    rho_vals = calc_rho_bar(k, P)

    return sigma_psi, q_vals, X_vals, Y_vals, eta_E, sigma_0, sigma_0_0, D_vals, F_vals, G_vals, H_vals, I_vals, rho_vals


def calc_covariances_func_one_loop(k, P, Q, R):

    # Function to calculate the covariances X, Y for LPT+1loop power spectrum

    sigma_psi = calc_sigma_psi(k, P)
    sigma_22 = calc_sigma_22(k, Q)
    sigma_13 = calc_sigma_13(k, R)
    q_vals, X_lin, Y_lin = calc_X_and_Y(k, P, sigma_psi)
    X22, Y22 = calc_X22_and_Y22(k, Q, sigma_22)
    X13, Y13 = calc_X13_and_Y13(k, R, sigma_13)

    X_loop = X_lin + X22 + 2.0*X13
    Y_loop = Y_lin + Y22 + 2.0*Y13

    return sigma_psi, sigma_22, sigma_13, q_vals, X_loop, Y_loop

# Function definitions


def dosph(n, x, f, a, tilt=1.5):

    # Function to calculate spherical Bessel integrals returns q values and integral values
    # Parameters: n = order of Bessel, x = variable of the function, f = function integrating over, a = power law

    func = renorm*np.power(x, a)*f(x)

    return sph(x, nu=n, q=tilt, lowring=True)(func, extrap=True)


def calc_sigma_psi(k, P):

    # Function to calculate the variance of the displacement field \sigma^2{\psi} = \frac{1}{6\pi^2}\int{dkP_L\left(k)}

    q, I0 = dosph(0, k, P, -2)
    sigma_vals = npi2*I0/6.0

    sigma_func = interp(q, sigma_vals)
    sigma_0 = sigma_func(min(q))

    return sigma_0


def calc_X_and_Y(k, P, sigma):

    # Function to calculate X and Y where X = \frac{1}{2\pi^2}\int{dk[\frac{2}{3}-\frac{j_1(kq)}{kq}]P_L} and Y = \frac{1}{2\pi^2}\int{dk[6\frac{j_1(kq)}{kq}-j_0(kq)]P_L} and sigma is \sigma^2_{\psi}

    q, I0 = dosph(0, k, P, -2)
    _, I2 = dosph(2, k, P, -2)

    X = 2.0*sigma - npi2*(1.0/3.0)*(I0+I2)
    Y = npi2*I2

    return q, X, Y


def calc_eta(k, P):

    # Function to calculate eta^2_E = \frac{1}{6\pi^2}\int{dk}k^2P_L

    q, I0 = dosph(0, k, P, 0)
    eta_vals = npi2*I0/6.0

    eta_func = interp(q, eta_vals)
    eta_0 = eta_func(min(q))

    return eta_0


def calc_sigma_0(k, P):

    # Function to calculate sigma^2_0=\frac{1}{\left(2\pi\right)^3}\int{dk}e^{-ik\cdotq}P_L

    q, I0 = dosph(0, k, P, 0)
    sigma_vals = 0.5*npi2*I0

    sigma_func = interp(q, sigma_vals)
    sigma_0 = sigma_func(min(q))

    return sigma_vals, sigma_0


def calc_D_and_F(k, P):

    # Function to calculate D=\frac{1}{2\pi^2}\int{dk kj_3(kq)P_L} and F=-\frac{1}{2\pi^2}\int{dk kj_2(kq)/kqP_L}

    q, I1 = dosph(1, k, P, -1)
    _, I3 = dosph(3, k, P, -1)

    F = -npi2*(I1+I3)/10.0
    D = 0.5*npi2*I3

    return q, D, F


def calc_G_H_and_I(k, P):

    # Function to calculate G=\frac{1}{2\pi^2}\int{dk k^2P_Lj_4(kq)}, H=-\frac{1}{14\pi^2}\int{dk k^2P_L(j_2+j_4)} and I=\frac{1}{210\pi^2}\int{dk k^2P_L(7j_0+10j_2+3j_4)}

    q, I0 = dosph(0, k, P, 0)
    _, I2 = dosph(2, k, P, 0)
    _, I4 = dosph(4, k, P, 0)

    G = 0.5*npi2*I4
    I = npi2*(7.0*I0+10.0*I2+3.0*I4)/210.0
    H = -npi2*(I2+I4)/14.0

    return q, G, H, I


def calc_rho_bar(k, P):

    q, I1 = dosph(1, k, P, -1)

    rho_bar = 0.5*npi2*I1

    return rho_bar


def calc_sigma_22(k, Q):

    # Function to calculate the variance of the displacement field \sigma^2{\psi} = \frac{1}{6\pi^2}\int{dkQ_1\left(k)}

    q, I0 = dosph(0, k, Q, -2)

    sigma_func = interp(q, npi2*I0/6.0)
    sigma_0 = sigma_func(min(q))

    return (9.0/98.0)*sigma_0


def calc_sigma_13(k, R):

    # Function to calculate the variance of the displacement field \sigma^2{\psi} = \frac{1}{6\pi^2}\int{dkR_1\left(k)}

    q, I0 = dosph(0, k, R, -2)

    sigma_func = interp(q, npi2*I0/6.0)
    sigma_0 = sigma_func(min(q))

    return (5.0/21.0)*sigma_0


def calc_X11_and_Y11(k, P, sigma):

    # Function to calculate X_11 and Y_11 where X_11 = \frac{1}{2\pi^2}\int{dk[\frac{2}{3}-\frac{j_1(kq)}{kq}]P_L} and Y_11 = \frac{1}{2\pi^2}\int{dk[6\frac{j_1(kq)}{kq}-j_0(kq)]P_L} and sigma is \sigma^2_{\psi}

    q, I0 = dosph(0, k, P, -2)
    _, I2 = dosph(2, k, P, -2)

    X11 = 2.0*sigma - npi2*(1.0/3.0)*(I0+I2)
    Y11 = npi2*I2

    return X11, Y11


def calc_X22_and_Y22(k, Q, sigma):

    # Function to calculate X_22 and Y_22 where X_22 = \frac{1}{2\pi^2}\int{dk[\frac{2}{3}-\frac{j_1(kq)}{kq}]Q_1} and Y_22 = \frac{1}{2\pi^2}\int{dk[6\frac{j_1(kq)}{kq}-j_0(kq)]Q_1} and sigma is \sigma^2_{22}

    q, I0 = dosph(0, k, Q, -2)
    _, I2 = dosph(2, k, Q, -2)

    X22 = 2.0*sigma-(9.0/98.0)*(npi2*(1.0/3.0)*(I0+I2))
    Y22 = (9.0/98.0)*npi2*I2

    return X22, Y22


def calc_X13_and_Y13(k, R, sigma):

    # Function to calculate X_13 and Y_13 where X_13 = \frac{1}{2\pi^2}\int{dk[\frac{2}{3}-\frac{j_1(kq)}{kq}]R_1} and Y_13 = \frac{1}{2\pi^2}\int{dk[6\frac{j_1(kq)}{kq}-j_0(kq)]R_1} and sigma is \sigma^2_{13}

    q, I0 = dosph(0, k, R, -2)
    _, I2 = dosph(2, k, R, -2)

    X13 = 2.0*sigma-(5.0/21.0)*(npi2*(1.0/3.0)*(I0+I2))
    Y13 = (5.0/21.0)*npi2*I2

    return X13, Y13
