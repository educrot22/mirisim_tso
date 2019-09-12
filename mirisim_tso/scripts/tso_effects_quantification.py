#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
------------------------------- Post-MIRISim -------------------------------
Post-treatment of MIRISim simulations for exoplanets time-serie observations
(Martin-Lagarde | CEA-Saclay | 2019)
"""

from pdb import set_trace as stop
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy.misc import factorial
#from scipy.stats import skewnorm


def __f_2exp(x, a, b, c, d, e): # OK
    result = a + b * np.exp(c*x) + d * np.exp(e*x)
    return result

def __f_exp(x, a, b, c, d, e): # OK
    result = a * np.exp(b * x) + c
    return result

def __f_log(x, a, b, c, d, e): # OK
    result = a + b*np.log(x) + c
    return result

def __linear(x, a, b, c, d, e): # OK
    result = a*x + b
    return result

def __rene(x, a, b, c, d, e): # OK
    result = (a*b*x+d)/(b +(x-c)**2)
    return result

def __polynomial(x, a, b, c, d, e): # OK
    result = a * x**2 + b * x + c
    #result = 2 * (M - D * log(M))
    return result



def print_formula(a, y, x, fit_form=None,):
    if fit_form=="log":
        print("Function | ", y, "(", x, ") = ", a[0]," + ", a[1], "x log(", x, ") +", a[2], "\n\n")
    elif fit_form=="exp":
        print("Function | ", y, "(", x, ") = ", a[0],"x exp(", a[1], "x ", x, ") +", a[2], "\n\n")
    elif fit_form=="2exp":
        print("Function | ", y, "(", x, ") = ", a[0]," + ", a[1],"x exp(", a[2], "x ", x, ") +", a[3],"x exp(", a[4], "x ", x, ")\n\n")
    elif fit_form=="lin":
        print("Function | ", y, "(", x, ") = ", a[0],"x ", x, " +", a[1], "\n\n")
    elif fit_form=="rene":
        print("Function | ", y, "(", x, ") = (", a[0]," x ", a[1]," x ", x," +", a[3],") / (", a[1]," + (", x," - ", a[2],")^2 )\n\n")
    elif fit_form=="poly":
        print("Function | ", y, "(", x, ") = ", a[0],"x ", x, "^2 +", a[1], "x ", x, " +",a[2],"\n\n")

    return

FORMS_AVAILABLE = {
            "log": __f_log,
            "exp": __f_exp,
            "2exp":__f_2exp,
            "lin": __linear,
            "rene":__rene,
            "poly":__polynomial,
        }
# > FITTING LAWS ON EVOLUTION TIME-CONSTANTS
# points_alpha = np.array([(10, 7000),(1800, 950), (2700, 800), (3600, 700), (5400, 600), (7200, 500), (36000, 5)]) # format (seconds, DN/s)   (8, 235000), (10, 7000),
# format           (alpha1,    alpha2,    DN/s)
# points = np.array([(208.70864, 3908.8192, 5000),
#                    (318.31608, 4299.5148, 4000),
#                    (544.46991, 3597.6787, 3000),
#                    (349.31304, 351.66609, 2000),
#                    (6.4332590, 1534.3659, 950),
#                    (59.537896, 1936.8444, 800),
#                    (90.853190, 2271.0439, 700),
#                    (109.24579, 2935.6206, 600),
#                    (125.62074, 3895.2158, 500),
#                    (142.39287, 5013.8033, 400),
#                    (145.25584, 5684.4315, 300),
#                    (171.44743, 4214.9832, 200),
#                    (164.52354, 1907.2530, 100),
#                    (190.29736, 1276.8339, 10)])




print("-------------------------------------------------------------------------------------------------------")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
print("                                        RESPONSE DRIFT                                                 ")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")


# format           (alpha1,    alpha2,    amplitude1, amplitude2, DN/s)
points = np.array([(263.72184, 1764.5250, -2.80276,   -0.660989,   11.73  ),
                   (212.48096, 2594.7549, -10.651909, -1.9204869,  86.5420),
                   (218.83157, 5167.1494, -9.7317493, -6.0712761,  206.844),
                   (180.19635, 6218.6697, -10.380856, -13.237013,  293.535),
                   (160.51597, 5095.6286, -8.8234513, -19.301959,  396.934),
                   (131.33408, 3909.9096, -11.360950, -22.762100,  504.068),
                   (97.374279, 2932.9071, -16.097594, -24.884621,  615.245),
                   (81.968308, 2272.5938, -18.055352, -27.562353,  712.214),
                   (53.726743, 1939.6679, -34.056237, -30.244137,  789.476), #])
                   (36.328816, 1633.2084, -104.83194, -32.746924,  884.091)])





                        # get x and y vectors
alpha1     = points[:,0]
alpha2     = points[:,1]
amplitude1 = points[:,2]
amplitude2 = points[:,3]
signal     = points[:,4]

# calculate polynomial
# z = np.polyfit(alpha, np.exp(signal), 2)
# f = np.poly1d(z)

print("\n • Alpha 1:")
fit_form="lin"
a_alpha1, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, alpha1)#, bounds=(0, np.inf))
print_formula(a_alpha1, "alpha1", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
alpha1_new = FORMS_AVAILABLE[fit_form](signal_new, a_alpha1[0], a_alpha1[1], a_alpha1[2], a_alpha1[3], a_alpha1[4])

print("\n • Alpha 2:")
fit_form="rene"
priors = (250., 200., 1., 6000., 1.0)
a_alpha2, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, alpha2, p0=priors)#, bounds=(0, np.inf))
print_formula(a_alpha2, "alpha2", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
# = (250., 1., 10000., a_alpha2[3], a_alpha2[4])
alpha2_new = FORMS_AVAILABLE[fit_form](signal_new, a_alpha2[0], a_alpha2[1], a_alpha2[2], a_alpha2[3], a_alpha2[4])


print("\n • Ampltude 1:")
fit_form="exp"
priors = (-1., 5e-3, -11., 1.0, 1.0)
a_amp1, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, amplitude1, p0=priors)#, bounds=(0, np.inf))
print_formula(a_amp1, "alpha1_amp", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
amplitude1_new = FORMS_AVAILABLE[fit_form](signal_new, a_amp1[0], a_amp1[1], a_amp1[2], a_amp1[3], a_amp1[4])

print("\n • Amplitude 2:")
fit_form="lin"
a_amp2, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, amplitude2)#, bounds=(0, np.inf))
print_formula(a_amp2, "alpha2_amp", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
amplitude2_new = FORMS_AVAILABLE[fit_form](signal_new, a_amp2[0], a_amp2[1], a_amp2[2], a_amp2[3], a_amp2[4])

# calculate new x's and y's
#first_point = FORMS_AVAILABLE[fit_form](8, a[0], a[1], a[2], a[3], a[4])



fig = plt.figure(figsize=(20,15))
fig.suptitle('FITTING FOR RESPONSE DRIFT TIME CONSTANTS', fontsize=16)
ax = fig.add_subplot(2, 2, 1)
plt.plot(signal, alpha1, 'o', signal_new, alpha1_new, color='xkcd:navy')
plt.title('Alpha 1')
#plt.legend()
# plt.xlim([0,500])
# plt.ylim([0,1500])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Time constant (s)')
# # -- # -- # -- # -- # -- # -- # -- #
ax = fig.add_subplot(2, 2, 2)
plt.plot(signal, alpha2, 'o',signal_new, alpha2_new, color='xkcd:raspberry')
plt.title('Alpha 2')
#plt.legend()
# plt.xlim([0,500])
# plt.ylim([0,1500])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Time constant (s)')
# -- # -- # -- # -- # -- # -- # -- #
ax = fig.add_subplot(2, 2, 3)
plt.plot(signal, amplitude1, 'o', signal_new, amplitude1_new, color='xkcd:navy')
plt.title('Amlplitude 1')
#plt.legend()
# plt.xlim([0,500])
plt.ylim([amplitude1.min(), amplitude1.max()])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Amplitude (DN/s)')

# # -- # -- # -- # -- # -- # -- # -- #
ax = fig.add_subplot(2, 2, 4)
plt.plot(signal, amplitude2, 'o', signal_new, amplitude2_new, color='xkcd:raspberry')
plt.title('Amlplitude 2')
#plt.legend()
# plt.xlim([0,500])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Amplitude (DN/s)')
#fig.tight_layout()


fig.savefig('alpha1_curvefits.png', dpi=300, transparent=True)









# format           (alpha3,    amplitude3,  DN/s)
points = np.array([(10327.691, -59.095313,  1995.68), #
                   (6133.6167, -74.831790,  2999.83), #
                   (3643.5560, -88.179656,  3997.67), #
                   (2283.5902, -93.478588,  5001.08)]) #


alpha3     = points[:,0]
amplitude3 = points[:,1]
signal     = points[:,2]


print("\n • Alpha 3:")
fit_form="exp"
prior = (4e6, -1e-2, 1000, 1.0, 1.0)
a_alpha3, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, alpha3, p0=prior)#, bounds=(0, np.inf))
print_formula(a_alpha3, "alpha3", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
beta1_new = FORMS_AVAILABLE[fit_form](signal_new, a_alpha3[0], a_alpha3[1], a_alpha3[2], a_alpha3[3], a_alpha3[4])


print("\n • Ampltude 3:")
fit_form="poly"
a_amp3, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, amplitude3)#, bounds=(0, np.inf))
print_formula(a_amp3, "alpha3_amp", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
amplitude1_new = FORMS_AVAILABLE[fit_form](signal_new, a_amp3[0], a_amp3[1], a_amp3[2], a_amp3[3], a_amp3[4])









fig = plt.figure(figsize=(20,15))
fig.suptitle('FITTING FOR RESPONSE DRIFT TIME CONSTANTS', fontsize=16)
ax = fig.add_subplot(2, 2, 1)
plt.plot(signal, alpha3, 'o', signal_new, alpha3_new, color='xkcd:navy')
plt.title('Alpha 1')
#plt.legend()
# plt.xlim([0,500])
# plt.ylim([0,1500])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Time constant (s)')
# -- # -- # -- # -- # -- # -- # -- #
ax = fig.add_subplot(2, 2, 3)
plt.plot(signal, amplitude3, 'o', signal_new, amplitude3_new, color='xkcd:navy')
plt.title('Amlplitude 1')
#plt.legend()
# plt.xlim([0,500])
plt.ylim([amplitude3.min(), amplitude3.max()])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Amplitude (DN/s)')


fig.savefig('alpha2_curvefits.png', dpi=300, transparent=True)












print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
print("                                        IDLE RECOVERY                                                  ")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")


# format           (beta1,    amplitude1,  DN/s     idle time (s))
points = np.array([(4222.2953, 0.49987286,  395.889, 2011.9), #
                   (1616.8302, 0.92332726,  501.280, 2011.9), #
                   (1231.9964, 1.58553411,  612.218, 2011.9), #
                   (1202.1463,  2.6876096,  708.617, 2011.9), #
                   (1167.8603,  3.8616116,  786.016, 2011.9), #
                   (1093.3717,  5.3563373,  880.532, 2011.9), #
                   (1074.6852,  6.5590930,  985.975, 2011.9), # dataset designe for exo-strategy
                   (1090.34,    10.314696,  1090.34, 2011.9), #
                   (1178.04,    13.050965,  1178.04, 2011.9)]) #
                   # (, , , , )])


# 15:21:49 + 387*12*2.775(s)
#
# 19:30:08
#
# 18:08:45 + 387*12*2.775(s)
#
# 21:46:42


                        # get x and y vectors
beta1     = points[:,0]
# beta2     = points[:,1]
amplitude1 = points[:,1]
# amplitude2 = points[:,3]
signal     = points[:,2]

# calculate polynomial
# z = np.polyfit(beta, np.exp(signal), 2)
# f = np.poly1d(z)

print("\n • beta 1:")
fit_form="exp"
prior = (4e6, -1e-2, 1000, 1.0, 1.0)
a_beta1, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, beta1, p0=prior)#, bounds=(0, np.inf))
print_formula(a_beta1, "beta", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
beta1_new = FORMS_AVAILABLE[fit_form](signal_new, a_beta1[0], a_beta1[1], a_beta1[2], a_beta1[3], a_beta1[4])


print("\n • Ampltude 1:")
fit_form="poly"
a_amp1, b = curve_fit(FORMS_AVAILABLE[fit_form], signal, amplitude1)#, bounds=(0, np.inf))
print_formula(a_amp1, "beta_amp", "signal", fit_form)
signal_new = np.linspace(signal[0], signal[-1], 200)
amplitude1_new = FORMS_AVAILABLE[fit_form](signal_new, a_amp1[0], a_amp1[1], a_amp1[2], a_amp1[3], a_amp1[4])




fig = plt.figure(figsize=(20,15))
fig.suptitle('FITTING FOR IDLE TIME TIME CONSTANTS', fontsize=16)
ax = fig.add_subplot(2, 1, 1)
plt.plot(signal, beta1, 'o', signal_new, beta1_new, color='xkcd:navy')
plt.title('beta 1')
plt.legend()
# plt.xlim([0,500])
plt.ylim([beta1.min(), beta1.max()])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Time constant (s)')
# -- # -- # -- # -- # -- # -- # -- #
ax = fig.add_subplot(2, 1, 2)
plt.plot(signal, amplitude1, 'o', signal_new, amplitude1_new, color='xkcd:navy')
plt.title('Amlplitude 1')
plt.legend()
# plt.xlim([0,500])
# plt.ylim([0,1500])
plt.xlabel('Signal (DN/s)')
plt.ylabel('Amplitude (DN/s)')
fig.savefig('beta_curvefits.png', dpi=300, transparent=True)





print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
print("                                       ANNEAL RECOVERY                                                 ")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")


fit_form="2exp"
a_anneal=(0.0, 11.600852, -1/197.97132, 0.86786327, -1/991.76323)
print_formula(a_anneal, "signal", "t_0", fit_form)
time = np.linspace(0, 6000, 200)
signal_anneal_new = FORMS_AVAILABLE[fit_form](time, a_anneal[0], a_anneal[1], a_anneal[2], a_anneal[3], a_anneal[4])
# -- # -- # -- # -- # -- # -- # -- #
fig = plt.figure()
plt.plot(time, signal_anneal_new, color='xkcd:navy')
plt.title('Anneal dark current decay')
plt.ylabel('Signal (DN/s)')
plt.xlabel('Time (s)')
fig.savefig('gamma_anneal_decay.png', dpi=300, transparent=True)

print("-------------------------------------------------------------------------------------------------------")
