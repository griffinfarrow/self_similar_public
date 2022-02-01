import numpy as np
import unnorm as un
from transport_coeff import *
import input as inp

no_variation_of_coul = inp.no_variation_of_coul #logical switch to decide if Coulomb log should be calculated self-consistently
l_unmagnetised = inp.l_unmagnetised

def n(y):
    theta, h, eps, q, v = y
    numerator = (1 + (1/un.beta_0) - h**2)
    denominator = theta
    if (l_unmagnetised):
        return ((1.0)/theta)
    else:
        return (numerator/denominator)

def dnde(x, y):
    theta, h, eps, q, v = y
    first = (-2*h/theta)*dhde(x, y)
    second = (-((1+(1/un.beta_0)-h**2)/(theta**2))*dtde(x, y))
    if (l_unmagnetised):
        return (-1.0/theta**2) * dtde(x, y)
    else:
        return (first+second) #(-1.0/theta**2)*dtde(x,y) #(first + second)

def dhde(x, y):
    theta, h, eps, q, v = y
    if (inp.nernst):
        prefactor = ((un.B_param * alpha(y)/(n(y) * tau_hat(y))) - (un.A_param * un.C_param * (beta(y)**2)/(n(y) * tau_hat(y) * gamma(y))))
        coeff = (1.0/prefactor)
        value = (eps - ((un.C_param * beta(y) * q)/(gamma(y) * n(y) * theta * tau_hat(y))))
    else:
        multiplier = 0.0
        prefactor = ((un.B_param * alpha(y)/(n(y) * tau_hat(y))) - multiplier * (un.A_param * un.C_param * (beta(y)**2)/(n(y) * tau_hat(y) * gamma(y))))
        coeff = (1.0/prefactor)
        value = (eps - multiplier * ((un.C_param * beta(y) * q)/(gamma(y) * n(y) * theta * tau_hat(y))))
    return (coeff * value)

def depsde(x, y):
    theta, h, eps, q, v = y
    prefactor = ((un.H0 * un.u0)/(un.c * un.E0))
    value = (v - (x/2.0)) * dhde(x, y) + h * dvde(x, y)
    # removing the frozen in flow contribution to the magnetic field
    # value = (-x/2.0) * dhde(x, y) + 0.0 * (v * dhde(x, y) + h * dvde(x, y))
    return (prefactor * value)

def dtde(x, y):
    theta, h, eps, q, v = y
    prefactor = (1.0/(gamma(y) * n(y) * theta * tau_hat(y)))
    first = q
    if (inp.ettings):
        second = -un.A_param * beta(y) * theta * dhde(x, y)
    else:
        second = 0.0
    return (prefactor * (first + second))

def dqde(x, y):
    theta, h, eps, q, v = y
    first = 3.0 * n(y) * (v - (x/2.0)) * dtde(x, y)
    second = 2.0 * n(y) * theta * dvde(x, y)
    if (inp.ohmic):
        third = -(eps/(np.sqrt(2*np.pi)))*dhde(x, y)
    else:
        third = 0.0
    return (first + second + third)

def dvde(x, y):
    theta, h, eps, q, v = y
    return ((1.0/n(y)) * ((x/2.0) - v) * dnde(x, y))

def mag(y):
    theta, h, eps, q, v = y
    return (un.mag_0 * abs(h) * tau_hat(y))

def delta(y):
    wt = mag(y)
    return ((wt**4) + delta_1 * (wt**2) + delta_0)

def alpha(y):
    wt = mag(y)
    numerator = alpha_1 * (wt**2) + alpha_0
    return (1.0 - (numerator/delta(y)))

def beta(y):
    wt = mag(y)
    numerator = wt * (beta_1 * (wt)**2 + beta_0)
    return numerator/delta(y)

def gamma(y):
    wete = mag(y) #electron magnetisation
    witi = (1/un.z)*np.sqrt(2*un.me/un.mi)*wete #ion magnetisation
    gamma_e = (gamma_1 * wete**2 + gamma_0)/(delta(y))
    gamma_i = (2.0 * witi**2 + 2.645)/(witi**4 + 2.70 * witi**2  + 0.677)
    return (gamma_e + np.sqrt(2*un.me/un.mi)*(1/un.z**3)*gamma_i)

def tau_hat(y):
    theta, h, eps, q, v = y
    numerator = theta**(1.5)
    denominator = lambdafunc(y) * n(y)
    return (numerator/denominator)

def lambdafunc(y):
    # function definition the Coulomb logarithm
    theta, h, eps, q, v = y
    if (no_variation_of_coul):
        return 1
    else:
        multiplier = 0.1
        if ( 1 + multiplier * (1/un.log_lambda_0)*(np.log(theta/np.sqrt(n(y)))) <= 0):
            raise Exception("non ideal plasma")
        return (1 + multiplier * (1/un.log_lambda_0)*np.log(theta/np.sqrt(n(y))))

def dyde(x, y):
    F = np.zeros(5)
    F[0] = dtde(x, y)
    F[1] = dhde(x, y)
    F[2] = depsde(x, y)
    F[3] = dqde(x, y)
    F[4] = dvde(x, y)
    return F
