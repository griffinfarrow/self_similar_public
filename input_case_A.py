# Script with all inputs for self-similar code

import numpy as np

################################################################################
## Normalisation constants
################################################################################

Z = 1
mi = 1.67e-24
N0 = 7.5e23
T0 = 1000
H0 = np.sqrt(16.0 * np.pi * N0 * T0 * 1.6e-12)
B0 = 1.0e7#0.2 * H0

################################################################################
## Boundary conditions (check r(u) in calculate.py)
################################################################################

h_negative_inf = B0/H0#0.19#24e6/H0#B0/H0
h_positive_inf = B0/H0#0.0168#24e6/H0#0.2#B0/H0
n_negative_inf = 1.0 #1.0 + h_positive_inf**2
n_positive_inf = 0.2
v_positive_inf = 0.0
v_negative_inf = 0.0

################################################################################
## Shooting optimisation settings
################################################################################

optimise = True
# estimate is [theta, h, eps, q, v]
estimate = [ 1.9767688623319044 , 0.0327020475311061 , 0.12216164629528062 , 1.777977030750518 , 0.35559540846434595 ]
#[ 1.986871513626507 , 0.03268829095517273 , 0.12408430987833456 , 1.8070159655311218 , 0.36140319358526946 ]
#[ 1.9809037854592466 , 0.03238173207281784 , 0.13224933872072756 , 1.90381120016875 , 0.3807622462442202 ]
#[1.978912039676056,0.032124216827193554,0.1427924441432785,2.0323147418216054,0.40646294664685323]
tol_value = 2.0e-6
optimise_f_stop = 20
optimise_r_stop = -20
damping = 1.0
h_val = 1.0e-8

################################################################################
## Values for full calculation
################################################################################
min_value = -20
max_value = 20

################################################################################
## Behaviour switches
################################################################################

separate = False
no_variation_of_coul = False # switch to control if Coulomb log should be calculated self-consistently
thomson = True # switch to control if we should allow the Thomson effect (heat flow with Nernst)
ohmic = False # switch to control if we should allow Ohmic heating
conduction = True # switch to control if we should have heat conduction
ettings = False # switch to control if we should switch off Ettingshausen effect
l_unmagnetised = True

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
save_data_path = "data/velikovich_attempt_differentcoullog/"
plot_real_values = False
plot = True
plotxmin = -20
plotxmax = 20
time = 1.0e-9 # physical time in s to plot the profiles for
output_as_Chimera = False # output data in a suitable way for Chimera to use
stylefile = "C:/Users/gf715/OneDrive - Imperial College London/style_files/aidan.mplstyle"
