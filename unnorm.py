import numpy as np
from tabulate import tabulate
import input as inp

z = inp.Z
mi = inp.mi #ion mass in grams
c = 2.9979e10 #speed of light
e = 4.8032e-10 #charge
hbar = 1.0546e-27 #reduced planck's constant
me = 9.1094e-28 #electron mass in grams
C_e = 1.6e-12 #conversion factor between eV and energy units

N0 = inp.N0
T0 = inp.T0 * C_e # converting eV to ergs
B0 = inp.B0
B0_Tesla = B0 * 1.0e-4
H0 = np.sqrt(8.0 * np.pi * (1 + (1/z)) * N0 * T0)

beta_0 = N0*T0*(1+1/z)/((B0**2)/(8*np.pi))
log_lambda_0 = np.log((T0/(e*hbar))*np.sqrt((3*me)/(np.pi*N0)))

####
# Alternative definition for Coulomb logarithm
#log_lambda_0 = 24 - np.log(np.sqrt(N0)/(T0/C_e))
#r = 0.3
#log_lambda_0 = np.log((T0/(e*hbar))*np.sqrt((3*me)/(np.pi*N0))) * (1-r) + (24 - np.log(np.sqrt(N0)/(T0/C_e))) * r
####

tau_e0 = 0.75*np.sqrt(me/(2.0*np.pi))*(T0**1.5)/(z*(e**4)*N0*log_lambda_0)
mag_0 = ((e * H0)/(me * c))*tau_e0
true_mag_0 = ((e * B0)/(me * c)) * tau_e0
Q0 = np.sqrt(((N0*T0)**2)*T0*tau_e0/me)
u0 = Q0/(N0*T0)
eta0 = (1.0/u0)
E0 = np.sqrt(8*np.pi) * Q0 / (c * H0)
A_param = c * T0 * H0 * (1/u0)/(4*np.pi*e*Q0)
B_param = c * me * H0 * (1/u0)/(4*np.pi*(e**2)*E0*N0*tau_e0)
C_param = T0*(1/u0)/(e*E0)

# some alternative normalisations for the new scheme to check they are equivalent


if __name__== "__main__":
    print("ratio =", B0/H0)
    print(tabulate([["N0 (/cm3)", N0], \
    ["T0 (eV)", T0/C_e], \
    ["B0 (G)", B0], \
    ["B0 (T)", B0_Tesla], \
    ["H0", H0], \
    ["beta0", beta_0], \
    ["log_lambda_0", log_lambda_0], \
    ["tau_e0 (s)", tau_e0],\
    ["mag_0", mag_0], \
    ["u0 (cm/s)", u0], \
    ["true magnetisation", true_mag_0], \
    ["A", A_param], \
    ["B", B_param], \
    ["C", C_param], \
    ["Q0", Q0], \
    ["E0", E0]], \
     headers = ['param', 'value']))
