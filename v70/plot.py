import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

DP = True
TP = True
# Funktion für Linfits
def linfit(x, a, b):
    return a*x + b

if DP == True:
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
    print(D0_err)


    ### Allgemeiner Plot 

    plt.errorbar(t, D0, yerr=D0_err, errorevery=(10), fmt=".", ecolor="pink", label="Messung")

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


'''
S von DP, Fit 1 =  1.39+/-0.14  l/s
m =  -0.0408+/-0.0007    b =  0.01+/-0.07

S von DP, Fit 2 =  0.49+/-0.05  l/s
m =  -0.01441+/-0.00034    b =  -4.16+/-0.08

S von DP, Fit 3 =  0.41+/-0.04  l/s
m =  -0.011943+/-0.000008    b =  -4.116+/-0.004

S von DP, Fit 4 =  0.203+/-0.020  l/s
m =  -0.005971+/-0.000004    b =  -6.174+/-0.006


Leck 0.5 mbar S =  1.02+/-0.15  l/s
m =  0.0150+/-0.0004    b =  1.63+/-0.04

Leck 10 mbar S =  1.60+/-0.17  l/s
m =  0.47184+/-0.00022    b =  14.866+/-0.016

Leck 50 mbar S =  1.64+/-0.17  l/s
m =  2.4190+/-0.0010    b =  59.78+/-0.07

Leck 1, 100 mbar S =  1.38+/-0.14  l/s
m =  4.0717+/-0.0017    b =  116.88+/-0.12

Leck 2, 100 mbar S =  1.41+/-0.15  l/s
m =  4.1531+/-0.0017    b =  118.92+/-0.13
'''
###################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
# 100 - 10^-8 mbar: 30%

pE = 1.35E-5 # Enddruck
VT = ufloat(33, 33*0.1) # Volumen Drehschieberpumpe

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak1_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3


### Allgemeiner Plot 

plt.errorbar(t, D2, yerr=D2_err, errorevery=(5), fmt=".", ecolor="pink", label="Messung")

# Zeiten für LinFit: 0 - 15, 15 - 50, 50 - Ende 
# 0 - 15
params, cov = curve_fit(linfit, t[0:15], np.log((D2[0:15]-pE)/(D2[0]-pE)), sigma=D2_err[0:15], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x015 = np.linspace(0, 20)
plt.errorbar(x015, D2[0]*np.exp(linfit(x015, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
print("S von TP, Fit 1 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 15 - 50
params, cov = curve_fit(linfit, t[15:50], np.log((D2[15:50]-pE)/(D2[0]-pE)), sigma=D2_err[15:50], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x1550 = np.linspace(10, 60)
plt.errorbar(x1550, D2[0]*np.exp(linfit(x1550, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
print("S von TP, Fit 2 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 50 - Ende
params, cov = curve_fit(linfit, t[50:], np.log((D2[50:]-pE)/(D2[0]-pE)), sigma=D2_err[50:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x50 = np.linspace(40, t[-1])
plt.errorbar(x50, D2[0]*np.exp(linfit(x50, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
print("S von TP, Fit 3 = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

plt.yscale("log")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Evakuierungskurve.pdf")
plt.figure()

# Fehler für Leckratenmessung
# 100 - 10^-8 mbar: 30%
# 
# # 1e-4 mbar
pG = unp.uarray(1E-4, (1E-4)*0.3)
t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_1e4_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D2_err = D2*0.3
plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VT/pG * linfit_params[0]
print("Leck 1E-4 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Leck_1e4.pdf")
plt.figure()

# # 2e-4 mbar
pG = unp.uarray(2E-4, (2E-4)*0.3)
t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_2e4_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D2_err = D2*0.3
plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VT/pG * linfit_params[0]
print("Leck 2E-4 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Leck_2e4.pdf")
plt.figure()

### Leitwerte
# t, D1, D2 = usecols(2, 4, 5)

# Fehler für Evakuierung:
# 100 - 10^-8 mbar: 30%

pE = 1.35E-5 # Enddruck
VT = ufloat(33, 33*0.1) # Volumen Drehschieberpumpe

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak1_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.errorbar(t, D2, yerr=D2_err, errorevery=(5), fmt=".", ecolor="pink", label="Messung", color="blue")

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak2_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.errorbar(t, D2, yerr=D2_err, errorevery=(5), fmt=".", ecolor="pink", label="Messung", color="red")

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak3_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.errorbar(t, D2, yerr=D2_err, errorevery=(5), fmt=".", ecolor="pink", label="Messung", color="green")

plt.yscale("log")
plt.savefig("TP_alle_evaks.pdf")
###
plt.figure()
t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak1_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.scatter(t, D2, s=1, label="Messung", color="blue")

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak2_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.scatter(t, D2, s=1, label="Messung", color="red")

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak3_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3

plt.scatter(t, D2, s=1, label="Messung", color="green")

plt.yscale("log")
plt.savefig("TP_alle_evaks_scatter.pdf")