import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

# Funktion für Linfits
def linfit(x, a, b):
    return a*x + b


#### Drehschieberpumpe

# t und D0 = usecols(2, 3)

# Fehler für Evakuierung
# 1000 - 10 mbar: 0.3%
# 10 - 2*10^-3:    10%
pE = 0.060 # Enddruck
VD = ufloat(34, 34*0.1) # Volumen Drehschieberpumpe

t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Evak_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D0_err = np.concatenate((D0[D0>=10]*0.003, D0[D0<10]*0.1)) # D0_err für Plots (nicht nützlich weil kaum erkennbar)


### Allgemeiner Plot 

plt.errorbar(t, D0, yerr=D0_err, errorevery=(10), ecolor="pink", label="Messung")

# Zeiten für LinFit: 0 - 150, 150 - 250, 250 - Ende 
# 0 - 150
params, cov = curve_fit(linfit, t[0:150], np.log((D0[0:150]-pE)/(D0[0]-pE)), sigma=D0_err[0:150], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VD
x0150 = np.linspace(0, 170)
plt.errorbar(x0150, D0[0]*np.exp(linfit(x0150, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
print("S von DP, Fit 1 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 150 - 250
params, cov = curve_fit(linfit, t[150:250], np.log((D0[150:250]-pE)/(D0[0]-pE)), sigma=D0_err[150:250], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VD
x150250 = np.linspace(130, 270)
plt.errorbar(x150250, D0[0]*np.exp(linfit(x150250, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
print("S von DP, Fit 2 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 250 - Ende
params, cov = curve_fit(linfit, t[250:], np.log((D0[250:]-pE)/(D0[0]-pE)), sigma=D0_err[250:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VD
x250end = np.linspace(230, 600)
plt.errorbar(x250end, D0[0]*np.exp(linfit(x250end, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
print("S von DP, Fit 3 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 250 - Ende mit künstlichem Fit
S = -linfit_params[0]*0.5*VD
plt.errorbar(x250end, D0[0]*np.exp(linfit(x250end, unp.nominal_values(linfit_params[0])*0.5, unp.nominal_values(linfit_params[1])*1.5)), color="black", label="Fit 4",)
print("S von DP, Fit 4 = ", S, " l/s")
print("m = ", linfit_params[0]*0.5, "   b = ", linfit_params[1]*1.5)
print(" ")

plt.yscale("log")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()

plt.savefig("DP_Evakuierungskurve.pdf")
plt.figure()
print(" ")

################################

# Fehler für Leckratenmessung
# 0.5 mbar:        10%
# Sonst:          0.3%

# 0.5 mbar
pG = unp.uarray(0.5, 0.05)
t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_05mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D0_err = D0*0.1

plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VD/pG * linfit_params[0]
print("Leck 0.5 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()
plt.savefig("DP_Leck_05mbar.pdf")
plt.figure()

# 10 mbar
pG = unp.uarray(10, 10*0.03)
t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_10mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D0_err = D0*0.003

plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VD/pG * linfit_params[0]
print("Leck 10 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()
plt.savefig("DP_Leck_10mbar.pdf")
plt.figure()

# 50 mbar
pG = unp.uarray(50, 50*0.03)
t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_50mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D0_err = D0*0.003

plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VD/pG * linfit_params[0]
print("Leck 50 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()
plt.savefig("DP_Leck_50mbar.pdf")
plt.figure()

# 100 mbar 1
pG = unp.uarray(100, 100*0.03)
t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck1_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D0_err = D0*0.003

plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VD/pG * linfit_params[0]
print("Leck 1, 100 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()
plt.savefig("DP_Leck1_100mbar.pdf")
plt.figure()

# 100 mbar 2
pG = unp.uarray(100, 100*0.03)
t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck2_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D0_err = D0*0.003

plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VD/pG * linfit_params[0]
print("Leck 2, 100 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
plt.legend()
plt.savefig("DP_Leck2_100mbar.pdf")
plt.figure()
plt.show()
plt.show()

###################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
# 100 - 10^-8 mbar: 30%

# Fehler für Leckratenmessung
# 100 - 10^-8 mbar: 30%

### Leitwerte
# t, D1, D2 = usecols(2, 4, 5)

# Fehler für Evakuierung:
# 100 - 10^-8 mbar: 30%

