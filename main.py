# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from decimal import Decimal
import pandas as pd
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from subprocess import call
import os
#from helpers import *
import helpers as hp
import sys
import re

font = {'family' : 'normal',
        'size'   :  12}
matplotlib.rc('font', **font)

# Custom code starts here


# some initial contants with their errors
T_bs = np.array([30.0, 50.0, 80.0, 110.0]) + 273.15 # °K
P_1 = hp.physical(200.0, 0.25*200.0, 3) # Ohm
R_1 = hp.physical(1200.0, 120.0, 4) # Ohm
R_2 = hp.physical(389*10**-3, 0.1*389*10**-3, 3) # Ohm
I_Shunt	= 20.0 #hp.physical(20.0, 0.1, 3) # A
V_Shunt = 50.0*10**-3 #hp.physical(50.0*10**-3, 0.1*50.0*10**-3, 2) # V
L_1 = 0.05 # m
L_2 = 0.05 # m
d_1 = 2.03*10**-3 # m
d_2 = 7.04*10**-3 # m
rho_1 = lambda T: np.polyval(np.polyfit(np.array([293.0, 353.0]), np.array([1.68, 2.06])*10**-8, 1), T) # linear polyfit
rho_2 = lambda T: 44.0*10**-8 # Ohm meters
k_1 = lambda T: np.polyval(np.polyfit([250.0, 400.0], [401.0, 391.0], 1), T) # W/mK
k_2 = lambda T: np.polyval(np.polyfit([275.0, 400.0], [21.9, 26.6], 1), T) # W/mK
F_1 = np.pi*(0.5*d_1)**2 # m^2
F_2 = np.pi*(0.5*d_2)**2 # m^2


# linear polyfit of the reference table of the type k thermocouple
C = np.arange(-10, 21, 1) # °C
V = np.array([-0.392, -0.353, -0.314, -0.275, -0.236, -0.197, -0.157, -0.118, -0.079, -0.039, 0.000, 0.039, 0.079, 0.119, 0.158, 0.198, 0.238, 0.277, 0.317, 0.357, 0.397, 0.437, 0.477, 0.517, 0.557, 0.597, 0.637, 0.677, 0.718, 0.758, 0.798]) * 10**-3 # V
coeffs1 = hp.phpolyfit(V, C, 1)
p1 = lambda x: np.polyval(coeffs1, x)

hp.replace("thermoFitline", hp.fmt_fit(coeffs1, None, 'V_T'))

# order of the colors in the plots
colors = {	"30":  "b",	# blue
			"50":  "g",	# green
			"80":  "r",	# red
			"110": "c"	# cyan
}


Pi12VTable, Pi12ITable = {}, {}
PiI = []
PiV = []
Pi = []

# create plots for all T's separately
for T_b in T_bs:

	Temp = str(int(T_b - 273.15)) # T in °C

	#if Temp == "80" or Temp == "50" or Temp == "110":
	# fetch the measured data from the xlsx file


	# the error of the multimeter is 1.2% + 1 digit
	deltaScale = 8.5*10**-6
	deltaV_S, deltaV_p, deltaI_T = "0.5%", "0.5%", "1.0%"
	deltaMMV = r"$\left ( \pm 0.5 \% + 1 \right) mV$"
	deltaMMA = r"$\left ( \pm 1.0 \% + 3 \right) \mu A$"
	deltaI_T_tot = r"$\left ( \pm 1.0 \% + 11.5 \right) \mu A$"
	hp.replace("deltaMMV", deltaMMV)
	hp.replace("deltaMMA", deltaMMA)
	hp.replace("deltaScale", deltaScale)
	hp.replace("deltaI_T_tot", deltaI_T_tot)
	V_S = hp.fetch2('data/example_data_' + Temp + '.xlsx', 'V_Shunt [mV]', deltaV_S)
	I_T = hp.fetch2('data/example_data_' + Temp + '.xlsx', 'I_T [microA]', deltaI_T)
	V_p = hp.fetch2('data/example_data_' + Temp + '.xlsx', 'V_p [mV]', deltaV_p)

	# add an error of 1 digit at the end
	V_S = hp.addError(V_S, 1.0)
	V_p = hp.addError(V_p, 1.0)
	I_T = hp.addError(I_T, 8.5 + 3.0)

	V_S = V_S*10**-3
	I_T = I_T*10**-6
	V_p = V_p*10**-3

	#else:
		# fetch the measured data from the csv file
		#V_S = hp.fetch('data/example_data_' + Temp + '.csv', 'V_Shunt [V]', 1.0*10**-3)
		#I_T = hp.fetch('data/example_data_' + Temp + '.csv', 'I_T [A]', deltaI_T)
		#V_p = hp.fetch('data/example_data_' + Temp + '.csv', 'V_p [V]', 1.0*10**-3)

	# do some calculation
	n = int(V_p.size/2) # n=4, half of the number of measurements
	V_T = R_2*I_T
	I = V_S*I_Shunt/V_Shunt
	V_p_mean = 0.5*(np.abs(np.flipud(V_p[0:n])) + np.abs(V_p[n:]))
	I_mean = 0.5*(np.abs(np.flipud(I[0:n])) + np.abs(I[n:]))
	deltaT = p1(V_T)
	deltaTminus = np.flipud(deltaT[0:n])
	deltaTplus = deltaT[n:]
	R_tot = L_1*rho_1(T_b)/F_1 + L_2*rho_2(T_b)/F_2
	Pi12V = 0.5*V_p_mean*(deltaTplus - deltaTminus)/(deltaTplus + deltaTminus)
	Pi12I = ((k_1(T_b)*F_1)/L_1 + (k_2(T_b)*F_2)/L_2)*(deltaTplus - deltaTminus)/(2*I_mean)
	PiI.append(np.mean(Pi12I))
	PiV.append(np.mean(Pi12V))
	Pi.append(np.mean(np.concatenate((Pi12V, Pi12I))))

	# some fitlines for the plots
	coeffs2 = hp.phpolyfit(I_mean, deltaTplus - deltaTminus, 1)
	p2 = lambda x: np.polyval(coeffs2, x)
	x2 = np.linspace(0, 18, 100)

	coeffs3 = hp.phpolyfit(I_mean, deltaTplus + deltaTminus, 2)
	p3 = lambda x: np.polyval(coeffs3, x)
	x3 = np.linspace(0, 18, 100)

	coeffs4 = hp.phpolyfit(I, deltaT, 2)
	p4 = lambda x: np.polyval(coeffs4, x)
	x4 = np.linspace(-18, 18, 100)

	### PLOTS ###

	# I - Pi12V
	plt.figure(0)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(Pi12V)*10**3, xerr=hp.stddev(I_mean), yerr=hp.stddev(Pi12V)*10**3, label=r'' + Temp + '°C')
	plt.ylabel(r'$\Pi_{12}^{(V)}$ [mV]')

	# I - Pi12I
	plt.figure(1)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(Pi12I)*10**3, xerr=hp.stddev(I_mean), yerr=hp.stddev(Pi12I)*10**3, label=r'' + Temp + '°C')
	plt.ylabel(r'$\Pi_{12}^{(I)}$ [mV]')

	# I - dT
	plt.figure(2)
	plt.errorbar(hp.nominal(I), hp.nominal(deltaT), xerr=hp.stddev(I), yerr=hp.stddev(deltaT), label=r'' + Temp + '°C', fmt='o', color=colors[Temp])
	plt.plot(x4, hp.nominal(p4(x4)), label=r'$f_{' + Temp + '^{\circ C}}(I)$')
	plt.ylabel(r'$\Delta T^{\pm}$ [°K]')
	hp.replace("dtpmfit"+Temp, hp.fmt_fit(hp.nominal(coeffs4), None, 'I'))

	# I - (dT+ - dT-)
	plt.figure(3)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(deltaTplus - deltaTminus), xerr=hp.stddev(I_mean), yerr=hp.stddev(deltaTplus - deltaTminus), label=r'' + Temp + '°C', fmt='o', color=colors[Temp])
	plt.plot(x2, hp.nominal(p2(x2)), label=r'$f_{' + Temp + '^{\circ C}}(I)$')
	plt.ylabel(r'$\Delta T^{+} - \Delta T^{-}$ [°K]')
	plt.xlim(left=0)
	plt.ylim(bottom=0)
	hp.replace("dtpmdtmfit"+Temp, hp.fmt_fit(hp.nominal(coeffs2), None, 'I'))

	# I - (dT+ + dT-)
	plt.figure(4)
	plt.errorbar(hp.nominal(I_mean), hp.nominal(deltaTplus + deltaTminus), xerr=hp.stddev(I_mean), yerr=hp.stddev(deltaTplus + deltaTminus), fmt='o', label=r'' + Temp + '°C', color=colors[Temp])
	plt.plot(x3, hp.nominal(p3(x3)), label=r'$f_{' + Temp + '^{\circ C}}(I)$')
	plt.ylabel(r'$\Delta T^{+} + \Delta T^{-}$ [°K]')
	hp.replace("dtppdtmfit"+Temp, hp.fmt_fit(hp.nominal(coeffs3), None, 'I'))

	# I - V_p
	plt.figure(5)
	plt.errorbar(hp.nominal(I), hp.nominal(V_p)*10**3, xerr=hp.stddev(I), yerr=hp.stddev(V_p)*10**3, label=r'' + Temp + '°C')
	plt.ylabel(r'$V_p$ [mV]')

	#plt.figure(6)
	#plt.errorbar(hp.nominal(T_bs), hp.nominal(Pi12I)*10**3, yerr=hp.stddev(Pi12I)*10**3, label=r'Pi12I at ' + Temp + 'C')
	#plt.errorbar(hp.nominal(T_bs), hp.nominal(Pi12V)*10**3, yerr=hp.stddev(Pi12V)*10**3, label=r'Pi12V at ' + Temp + 'C')

	# attributes of all plots
	for i in [0, 1, 2, 3, 4, 5]:
		plt.figure(i)
		plt.legend(loc="best")
		plt.grid(True)
		plt.xlabel(r'$I$ [A]')

	V_S = np.concatenate((np.array([r"$V_S [mV]$"]), V_S*10**3))
	I_T = np.concatenate((np.array([r"$I_T [\mu A]$"]), I_T*10**6))
	V_p = np.concatenate((np.array([r"$V_p [mV]$"]), V_p*10**3))
	I   = np.concatenate((np.array([r"$I [A]$"]), I))
	dT  = np.concatenate((np.array([r"$\Delta T^{\pm} [^{\circ}K]$"]), deltaT))
	arr = np.array([I, V_S, I_T, V_p, dT]).T
	hp.replace("table"+Temp, arr)

	Pi12VTable[Temp] = np.concatenate((np.array([r"$T_b=" + Temp + r"^{\circ}C$", r"$\Pi_{12}^{(V)}\, [mV]$"]), Pi12V*10**3))
	Pi12ITable[Temp] = np.concatenate((np.array([r"$T_b=" + Temp + r"^{\circ}C$", r"$\Pi_{12}^{(I)}\, [mV]$"]), Pi12I*10**3))

#T_bs -= 273.15
T_bs = np.array([
	hp.physical(30, 1, 2),
	hp.physical(50, 1, 2),
	hp.physical(80, 1, 2),
	hp.physical(110, 1, 3)
	])

# T - Pi
# Literature values
# https://www.researchgate.net/profile/Fabrice_Clerot/post/What_is_Seebeck_coefficient_of_copper_and_constantan_alloy_of_45_nickel_and_55Cu_at_room_temperature/attachment/59d62c0dc49f478072e9dcb8/AS%3A273540627009560%401442228576018/download/Guan+-+Seebeck+measurement+apparatus.pdf
S = lambda x: 4.37184 + 0.1676*x - 1.84371*10**(-4)*x**2 + 1.2244*10**(-7)*x**3 - 4.47618*10**(-11)*x**4
P = lambda x: x*S(x)
TT = np.linspace(300, 390, 100)
plt.figure(6)
#T_bs += 273.15
plt.errorbar(hp.nominal(T_bs), hp.nominal(PiI)*10**3, xerr=hp.stddev(T_bs), yerr=hp.stddev(PiI)*10**3, label=r'$\Pi_{12}^{(I)}$')
plt.errorbar(hp.nominal(T_bs), hp.nominal(PiV)*10**3, xerr=hp.stddev(T_bs), yerr=hp.stddev(PiV)*10**3, label=r'$\Pi_{12}^{(V)}$')
plt.errorbar(hp.nominal(T_bs), hp.nominal(Pi)*10**3, xerr=hp.stddev(T_bs), yerr=hp.stddev(Pi)*10**3, label=r'$\Pi_{12}$')
plt.plot(TT - 273.15, P(TT)*10**-3, label=r'Literature value: $\Pi(T)$')
plt.xlabel(r'$T \, [^{\circ}C]$')
plt.ylabel(r'$\Pi \, [mV]$')
plt.legend(loc="best")
#T_bs -= 273.15

T_bs = np.concatenate((np.array([r"$T_b \, [^{\circ}C]$"]), T_bs))
PiI = np.concatenate((np.array([r"$<\Pi_{12}^{(I)}>_I \, [mV]$"]), np.array(PiI)*10**3))
PiV = np.concatenate((np.array([r"$<\Pi_{12}^{(V)}>_I \, [mV]$"]), np.array(PiV)*10**3))
Pi = np.concatenate((np.array([r"$<\Pi_{12}>_I \, [mV]$"]), np.array(Pi)*10**3))
arr = np.array([T_bs, PiI, PiV, Pi]).T
hp.replace("PeltierTable", arr)

I = np.concatenate((np.array([r"", r"$I [A]$"]), np.array([4, 8, 12, 16])))
arr = np.array([I,
	Pi12VTable["30"],
	Pi12VTable["50"],
	Pi12VTable["80"],
	Pi12VTable["110"],
]).T
hp.replace("Pi12VTable", arr)

arr = np.array([I,
	Pi12ITable["30"],
	Pi12ITable["50"],
	Pi12ITable["80"],
	Pi12ITable["110"],
]).T
hp.replace("Pi12ITable", arr)

# save the plots im format .fmt
fmt = "eps"
plt.figure(0)
plt.savefig("plots/Pi12V_vs_I." + fmt)

plt.figure(1)
plt.savefig("plots/Pi12I_vs_I." + fmt)

plt.figure(2)
plt.savefig("plots/dT_vs_I." + fmt)

plt.figure(3)
plt.savefig("plots/dTp-dTm_vs_I." + fmt)

plt.figure(4)
plt.savefig("plots/dTp+dTm_vs_I." + fmt)

plt.figure(5)
plt.savefig("plots/Vp_vs_I." + fmt)

plt.figure(6)
plt.savefig("plots/Pi_vs_T." + fmt)

# replace the markers in main.tex with the formatted variables/fitlines/tables etc...
hp.replace("variable1", hp.physical(12345678923.3, 0.00000345678923, 3))
hp.replace("variable2", hp.physical(69.911212, 6.721212, 4))
hp.replace("variable3", hp.physical(0.00992123123, 0.00095123123, 3))
hp.replace("variable4", Pi12I[0])
hp.replace("variable5", Pi12V[0])
hp.replace("variable6", 0.0)
hp.replace("variable7", 0)
hp.replace("variable8", hp.physical(20.0, 0.000001, 3))
hp.replace("variable9", hp.physical(20.0, 0.000003, 3))
hp.replace("variable10", 10**23)
hp.replace("polyfit1", "f(x) = " + hp.fmt_fit(coeffs1, 1))
hp.replace("polyfit2", hp.fmt_fit(coeffs1, 2))
hp.replace("polyfit3", hp.fmt_fit([2, 5], 1))
hp.replace("polyfit4", hp.fmt_fit(coeffs1))


hp.replace("Name", "Roman Gruber")
hp.replace("Experiment", "Peltier Effect")
hp.replace("I_Shunt", hp.fmt_number(I_Shunt, 2))
hp.replace("V_Shunt", hp.fmt_number(V_Shunt*10**3, 2))
hp.replace("R_2", R_2*10**3)

# alternative measurement for T_b = 100°C
deltaV_S, deltaV_p, deltaI_T = "0.5%", "0.5%", "1.0%"
V_S = hp.fetch2('data/example_data_110_old.xlsx', 'V_Shunt [mV]', deltaV_S)
I_T = hp.fetch2('data/example_data_110_old.xlsx', 'I_T [microA]', deltaI_T)
V_p = hp.fetch2('data/example_data_110_old.xlsx', 'V_p [mV]', deltaV_p)
V_S = hp.addError(V_S, 1.0)
V_p = hp.addError(V_p, 1.0)
I_T = hp.addError(I_T, 8.5 + 3.0)
V_S = V_S*10**-3
I_T = I_T*10**-6
V_p = V_p*10**-3
I = V_S*I_Shunt/V_Shunt
V_S = np.concatenate((np.array([r"$V_S [mV]$"]), V_S*10**3))
I_T = np.concatenate((np.array([r"$I_T [\mu A]$"]), I_T*10**6))
V_p = np.concatenate((np.array([r"$V_p [mV]$"]), V_p*10**3))
I   = np.concatenate((np.array([r"$I [A]$"]), I))
arr = np.array([I, V_S, I_T, V_p]).T
hp.replace("table1102", arr)

# alternative measurement for T_b = 100°C
deltaV_S, deltaV_p, deltaI_T = "0.5%", "0.5%", "1.0%"
V_S = hp.fetch2('data/example_data_80_old.xlsx', 'V_Shunt [mV]', deltaV_S)
I_T = hp.fetch2('data/example_data_80_old.xlsx', 'I_T [microA]', deltaI_T)
V_p = hp.fetch2('data/example_data_80_old.xlsx', 'V_p [mV]', deltaV_p)
V_S = hp.addError(V_S, 1.0)
V_p = hp.addError(V_p, 1.0)
I_T = hp.addError(I_T, 8.5 + 3.0)
V_S = V_S*10**-3
I_T = I_T*10**-6
V_p = V_p*10**-3
I = V_S*I_Shunt/V_Shunt
V_S = np.concatenate((np.array([r"$V_S [mV]$"]), V_S*10**3))
I_T = np.concatenate((np.array([r"$I_T [\mu A]$"]), I_T*10**6))
V_p = np.concatenate((np.array([r"$V_p [mV]$"]), V_p*10**3))
I   = np.concatenate((np.array([r"$I [A]$"]), I))
arr = np.array([I, V_S, I_T, V_p]).T
hp.replace("table802", arr)

hp.compile()

#plt.show()