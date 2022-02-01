import numpy as np
import equations as eq
from newton import newtonRaphson2
import scipy.integrate as spi
from tabulate import tabulate
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interpol
import unnorm as un
import scipy.optimize as spo
import input as inp

mpl.rcParams['axes.formatter.useoffset'] = False
#plt.style.use(inp.stylefile)

def stick_together(x_f, x_b, forward, backward):
    x = np.concatenate(((np.flip(x_b)[:-1], x_f)))
    theta = np.concatenate((np.flip(backward[0,:])[:-1], forward[0,:]))
    h = np.concatenate((np.flip(backward[1,:])[:-1], forward[1,:]))
    eps = np.concatenate((np.flip(backward[2,:])[:-1], forward[2,:]))
    q = np.concatenate((np.flip(backward[3,:])[:-1], forward[3,:]))
    v = np.concatenate((np.flip(backward[4,:])[:-1], forward[4,:]))
    y = np.array([theta, h, eps, q, v])
    n = eq.n(y)
    return x, theta, h, eps, q, v, n

def r(val):
    r = np.zeros(np.size(val))
    start = 0.0
    stop_r = inp.optimise_r_stop
    stop_f = inp.optimise_f_stop
    if (inp.high_resolution):
        solution_r = spi.solve_ivp(eq.dyde, [start, stop_r], val, method = inp.method, max_step = inp.high_res_value)
        solution_f = spi.solve_ivp(eq.dyde, [start, stop_f], val, method = inp.method, max_step = inp.high_res_value)
    else:
        solution_r = spi.solve_ivp(eq.dyde, [start, stop_r], val, method = inp.method)
        solution_f = spi.solve_ivp(eq.dyde, [start, stop_f], val, method = inp.method)
    Y_r = solution_r.y
    Y_f = solution_f.y
    y_r = Y_r[:,-1]
    y_f = Y_f[:,-1]
    h_r = y_r[1]
    h_f = y_f[1]
    n_r = eq.n(y_r)
    n_f = eq.n(y_f)
    r[0] = y_r[1] - inp.h_negative_inf
    r[1] = n_r - inp.n_negative_inf
    #r[1] = y_r[0] - inp.T_negative_inf
    r[2] = y_f[1] - inp.h_positive_inf
    r[3] = n_f - inp.n_positive_inf
    #r[3] = y_f[0] - inp.T_positive_inf
    r[4] = y_r[4] - inp.v_positive_inf
    return r

def calc_deriv(x, f):
    deriv = np.zeros_like(f)
    for i in range(np.size(x)-1):
        deriv[i] = (f[i+1] - f[i])/(x[i+1] - x[i])
    return deriv

print("-------------------------------------------------------------------------")
print("-------------------------------------------------------------------------")
print("using a density jump of ", inp.n_negative_inf)
print("\nusing a field of ", inp.h_positive_inf)
if (inp.high_resolution):
    print("\nusing a (high) resolution of ", inp.high_res_value)

print("-------------------------------------------------------------------------")
print("-------------------------------------------------------------------------")
print("NORMALISATION VALUES")
print(tabulate([["N0 (/cm3)", un.N0], \
["T0 (eV)", un.T0/un.C_e], \
["B0 (G)", un.B0], \
["B0 (T)", un.B0_Tesla], \
["beta_0", un.beta_0], \
["H0", un.H0]], headers = ['param', 'value']))

print("\nNOW BEGINNING SHOOTING PROBLEM")

print("\nMAGNETISED")

u = np.array(inp.estimate)
if (inp.optimise):
    u = newtonRaphson2(r, u, tol = inp.tol_value, full_output = True)

print("\noptimized u ", "[", u[0], ",", u[1], ",", u[2], ",", u[3], ",", u[4], "]")
if (not inp.high_resolution):
    print("\noptimized r ",  r(u))
print("\nNOW CALCULATING FULL SOLUTION")

xStart = 0.0
xStop = inp.min_value
if (inp.high_resolution):
    print("calculating in high resolution")
    solution_r = spi.solve_ivp(eq.dyde, [xStart, xStop], u, method = inp.method, max_step = inp.high_res_value)
else:
    print("calculating in low resolution")
    solution_r = spi.solve_ivp(eq.dyde, [xStart, xStop], u, method = inp.method)
x_r = solution_r.t
Y_r = solution_r.y

xStop = inp.max_value
if (inp.high_resolution):
    solution_f = spi.solve_ivp(eq.dyde, [xStart, xStop], u, method = inp.method, max_step = inp.high_res_value)
else:
    solution_f = spi.solve_ivp(eq.dyde, [xStart, xStop], u, method = inp.method)
x_f = solution_f.t
Y_f = solution_f.y


x, theta, h, eps, q, v, n = stick_together(x_f, x_r, Y_f, Y_r)

print(tabulate([["h", h[0], h[-1]],\
 ["theta",theta[0], theta[-1]], \
["q", q[0], q[-1]],\
 ["n", n[0], n[-1]], \
 ["v", v[0], v[-1]], \
 ["eps", eps[0], eps[-1]]],\
 headers = ['param', '- inf', '+ inf']))

if (inp.plot):
    fig, ax = plt.subplots(3, 1, figsize = (12, 6))
    ax[0].plot(x, h, 'b-')
    ax[0].set_ylabel(r"$h$", color = 'b')
    ax[0].set_xlabel(r"$\eta$")
    ax1 = ax[0].twinx()
    ax1.set_xlim([inp.plotxmin,inp.plotxmax])
    ax1.plot(x, n, 'r--')
    ax1.set_ylabel(r"$n$", color = 'r')
    ax[1].plot(x, theta)
    ax[1].set_xlabel(r"$\eta$")
    ax[1].set_ylabel(r"$\theta$")
    ax[1].set_xlim([inp.plotxmin,inp.plotxmax])
    ax[2].plot(x, v)
    ax[2].set_xlabel(r"$\eta$")
    ax[2].set_ylabel(r"$v$")
    ax[2].set_xlim([inp.plotxmin,inp.plotxmax])
    fig.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "all.png")

if (inp.plot_real_values):
    # plots the self-similar code unnormalised to its real values
    x_cm = x * np.sqrt(inp.time) / un.eta0
    x_m = x_cm * 0.01
    x_mm = x_m * 1000.0
    print("first xcm = ", x_cm[0])
    print("final xcm = ", x_cm[-1])
    num_dens = n * un.N0
    mass_dens_kgm3 = num_dens * 1.0e6 * 1.67e-27
    T_erg = theta * un.T0
    T_eV = T_erg/un.C_e
    v_cms = un.u0 * v / np.sqrt(inp.time)
    v_ms = v_cms * 0.01
    B_G = h * un.H0

    plt.figure()
    plt.plot(x_m, num_dens)
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$n$")
    plt.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "n.png")

    plt.figure()
    plt.plot(x_m, mass_dens_kgm3)
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$\rho$")
    plt.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "rho.png")

    plt.figure()
    plt.plot(x_m, v_ms)
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$v$ (m/s)")
    plt.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "v.png")

    plt.figure()
    plt.plot(x_m, T_eV)
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$T$ (eV)")
    plt.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "T.png")

    plt.figure()
    plt.plot(x_mm, B_G/B_G[-1])
    plt.xlabel(r"$x$ (mm)")
    plt.ylabel(r"$B/B_0$")
    plt.xlim(-0.1, 0.15)
    plt.tight_layout()
    if (inp.save_data):
        plt.savefig(inp.save_data_path + "B.png")

if (inp.save_data):
    print("\nNOW SAVING DATA")
    path_for_data = inp.save_data_path
    if (inp.high_resolution):
        path_for_data = path_for_data + "high_res/"
    np.savetxt(path_for_data + "x.txt", x)
    np.savetxt(path_for_data + "theta.txt", theta)
    np.savetxt(path_for_data + "n.txt", n)
    np.savetxt(path_for_data + "v.txt", v)
    np.savetxt(path_for_data + "eps.txt", eps)
    np.savetxt(path_for_data + "h.txt", h)
    np.savetxt(path_for_data + "q.txt", q)

# if (inp.check_energy_conservation):
#     print("\nCHECKING ENERGY CONSERVATION")
#     x_cm = x * np.sqrt(inp.time) / un.eta0
#     x_m = x_cm * 0.01
#     num_dens = n * un.N0
#     n_m3 = num_dens * 1.0e6
#     T_erg = theta * un.T0
#     T_eV = T_erg/un.C_e
#     v_cms = un.u0 * v / np.sqrt(inp.time)
#     v_ms = v_cms * 0.01
#     B_T = h * un.H0 * 1.0e-4
#     mass_dens = num_dens * un.mi * 1000.0
#
#     # set up and plot initial values
#     example_T = np.zeros(np.size(x_m))
#     example_n = np.zeros(np.size(x_m))
#     example_B = np.zeros(np.size(x_m))
#     for i in range(np.size(x_m)):
#         if x_m[i] <= 0:
#             example_n[i] = 2.0 * un.N0 * 1.0e6
#             example_T[i] = un.T0 / un.C_e
#             example_B[i] = 0.0
#         else:
#             example_n[i] = un.N0 * 1.0e6
#             example_T[i] = un.T0 / un.C_e
#             example_B[i] = un.H0 * 1.0e-4
#     example_energy = (((example_B**2)/(2*4.0*np.pi*1.0e-7)) + (1.5 * example_n * 1.6e-19 * example_T))
#
#     plt.figure(figsize = (12, 6))
#     plt.plot(x_m, (((B_T**2)/(2*4.0*np.pi*1.0e-7)) + (1.5 * n_m3 * 1.6e-19 * T_eV)), label = "late time")
#     plt.plot(x_m, example_energy, label = "initially")
#     plt.ylabel(r"$E$")
#     plt.xlabel(r"$x$")
#     plt.legend()
#     plt.tight_layout()

if (inp.output_as_Chimera):
    print("=======================================")
    print("Parameters for Chimera")
    print("=======================================")
    # convert to a suitable output datatype for Chimera
    # Chimera temperatures are in eV
    x_cm = x * np.sqrt(inp.time) / un.eta0
    x_m = x_cm * 0.01
    num_dens = n * un.N0
    T_erg = theta * un.T0
    T_eV = T_erg/un.C_e
    v_cms = un.u0 * v / np.sqrt(inp.time)
    v_ms = v_cms * 0.01
    B_T = h * un.H0 * 1.0e-4
    mass_dens = num_dens * un.mi * 1000.0

    print("T_HS (eV) = ", T_eV[0])
    print("rho_HS (kg/m3) = ", mass_dens[0])
    print("B_init (T) = ", B_T[-1])
    print("v_LHS (m/s) = ", v_ms[0])
    print("v_RHS (m/s) = ", v_ms[-1])
    print("peak velocity (m/s) = ", np.max(abs(v_ms)))
print("-------------------------------------------------------------------------")
print("-------------------------------------------------------------------------")

plt.show()
