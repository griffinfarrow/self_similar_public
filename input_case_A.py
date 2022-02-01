# Script with all inputs for self-similar code

import numpy as np

################################################################################
## Normalisation constants
################################################################################
"""
This section specifies the physical values to the use for e.g. density,
temperature and field. These are the values that things are normalised to

Z - The atomic number of the plasma. Code is only tested for Z = 1, but nothing
in principle wouldn't work for Z != 1
mi - The mass of the ions in the plasma
N0 - the reference density in /cm3
T0 - the reference temperature in eV
H0 - the reference value of magnetic field. DO NOT SET MANUALLY. Calculated from
     N0 and T0
B0 - the physical value of magnetic field (in G)
"""

Z = 1
mi = 1.67e-24
N0 = 7.5e23
T0 = 1000
H0 = np.sqrt(16.0 * np.pi * N0 * T0 * 1.6e-12)
B0 = 1.0e7#0.2 * H0

################################################################################
## Boundary conditions (check r(u) in calculate.py)
################################################################################

"""
Boundary conditions that we are trying to satisfy
The code will not necessarily satisfy all of these (can only satisfy 5 of them)
To decide which are satisfied, look at the function r(u) in calculate.py. Editing
that function edits the boundary conditions that are actually enforced
CAUTION: These boundary conditions are applied to the normalised values

h_negative_inf - the value of h at eta = -inf
h_positive_inf - the value of h at eta = inf
n_negative_inf - the value of n at eta = -inf
n_positive_inf - the value of n at eta = inf
v_positive_inf - the value of v at eta = inf
v_negative_inf - the value of v at eta = -inf
"""

h_negative_inf = B0/H0#0.19#24e6/H0#B0/H0
h_positive_inf = B0/H0#0.0168#24e6/H0#0.2#B0/H0
n_negative_inf = 1.0 #1.0 + h_positive_inf**2
n_positive_inf = 0.2
v_positive_inf = 0.0
v_negative_inf = 0.0

################################################################################
## Shooting optimisation settings
################################################################################

"""
Defining parameters for the shooting problem opitmisation
optimise - logical switch. If False, then no shooting problem optimisation is
           actually carried out
estimate - The estimate for the values on the eta = 0 interface. A poor initial
           guess will lead to the answer not converging. A list of reals
tol_value - The tolerance at which we should stop attempting shooting optimisation
            Basically if |r(latest_guess) - r(previous_guess)| < tol_value, will
            terminate
optimise_f_stop - The value of eta at which to stop optimisation on the forward
                  calculation (to eta = inf)
optimise_r_stop - The value of eta at which to stop optimisation on the backward
                  calculation (to eta = -inf)
damping - A float, slows down the Newton-Raphson calculation. Basically becomes
          x_{n+1} = x_n - damping*f(x_n)/f'(x_n) instead. Can help to find
          troublesome roots
h_val - The size of the step used in estimating the Jacobian by finite differences
        Edit with caution
"""

optimise = True
# estimate is [theta, h, eps, q, v]
estimate = [ 1.9767688623319044 , 0.0327020475311061 , 0.12216164629528062 , 1.777977030750518 , 0.35559540846434595 ]
tol_value = 2.0e-6
optimise_f_stop = 20
optimise_r_stop = -20
damping = 1.0
h_val = 1.0e-8

################################################################################
## Values for full calculation
################################################################################

"""
The maximum and minimum values of eta for the calculation of the solution
"""

min_value = -20
max_value = 20

################################################################################
## Behaviour switches
################################################################################

"""
Logical Behaviour to control the behaviour of the system
All are pretty clear, except for

no_variation_of_coul - for historical reasons, this is basically a double negative.
                       If no_variation_of_coul = False, then the Coulomb logarithm
                       is allowed to vary with eta.

l_unmagnetised - v. rarely used and will often break the code. Changes the pressure
                 balance to just use thermal pressure rather than magnetic
                 pressure also
"""

no_variation_of_coul = False # switch to control if Coulomb log should be calculated self-consistently
thomson = True # switch to control if we should allow the Thomson effect (heat flow with Nernst)
ohmic = False # switch to control if we should allow Ohmic heating
conduction = True # switch to control if we should have heat conduction
ettings = False # switch to control if we should switch off Ettingshausen effect
l_unmagnetised = True

################################################################################
## General switches
################################################################################

"""
Further behaviour switches
method - The method used for the ODE part of the shooting method. See
https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
for the various options

high_resolution - Switch to control whether we should limit the step size.
                  Rarely needed unless you're going to be calculating derivatives
                  of the solution (in which case it is helpful). CAUTION - very
                  slow

high_res_value - The value for the maximum step size in eta
"""

method = "BDF" # switch for integration method, "BDF" works best
high_resolution = False # calculate using a fixed small step size
high_res_value = 0.01

################################################################################
## Data output
################################################################################

"""
Controlling the data output
save_data - Switch to determine if data should be saved
save_data_path - Path to ALREADY EXISTING directory to save data into
plot_real_values - Logical switch. Whether to calculate the parameters in real
                   space (i.e. x and t rather than eta)
plot - Logical switch to determine if anything should be plot
plotxmin - minimum x value for the plot
plotxmax - Maximum x value for the plot
time - Real. The time at which to calculate the real values
output_as_Chimera - logical switch. Whether to produce the outputs in a suitable
                    way to be an initial condition for Chimera
stylefile - Path to a matplotlib stylefile for the plots
"""

save_data = True
save_data_path = "data/velikovich_attempt_differentcoullog/"
plot_real_values = False
plot = True
plotxmin = -20
plotxmax = 20
time = 1.0e-9 # physical time in s to plot the profiles for
output_as_Chimera = False # output data in a suitable way for Chimera to use
stylefile = "C:/Users/gf715/OneDrive - Imperial College London/style_files/aidan.mplstyle"
