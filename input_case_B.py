# Script with all inputs for self-similar code

import numpy as np

################################################################################
## Normalisation constants
################################################################################

Z = 1
mi = 1.67e-24
N0 = 1.0e20#5.0e19
T0 = 250
H0 = np.sqrt(16.0 * np.pi * N0 * T0 * 1.6e-12)
B0 = H0

################################################################################
## Boundary conditions (check r(u) in calculate.py)
################################################################################

h_negative_inf = 0.0
h_positive_inf = B0/H0
n_negative_inf = 5.0 #1.0 + h_positive_inf**2
n_positive_inf = 1.0
v_positive_inf = 0.0
v_negative_inf = 0.0

################################################################################
## Shooting optimisation settings
################################################################################

optimise = True
# estimate is [theta, h, eps, q, v]
estimate = [ 0.48131497574599347 , 0.5485494245422294 , 0.2670082679881751 , 0.14010830583385456 , 0.017432267790572605 ]
tol_value = 1e-5
optimise_f_stop = 30
optimise_r_stop = -30
damping = 1.0
h_val = 1.0e-4

################################################################################
## Values for full calculation
################################################################################
min_value = -30
max_value = 30

################################################################################
## Behaviour switches
################################################################################

no_variation_of_coul = False # switch to control if Coulomb log should be calculated self-consistently
ohmic = True # switch to control if we should allow Ohmic heating
ettings = True # switch to control if we should switch off Ettingshausen effect
l_unmagnetised = False
nernst = False

################################################################################
## General switches
################################################################################

method = "BDF" # switch for integration method, "BDF" works best
high_resolution = False # calculate using a fixed small step size
high_res_value = 0.01

################################################################################
## Data output
################################################################################
save_data = True
save_data_path = "data/unity_beta_temp_gradient/no_nernst/"
plot_real_values = True
plot = True
plotxmin = -30
plotxmax = 30
time = 1.0e-9 # physical time in s to plot the profiles for
output_as_Chimera = True # output data in a suitable way for Chimera to use
stylefile = "C:/Users/gf715/OneDrive - Imperial College London/style_files/aidan.mplstyle"
