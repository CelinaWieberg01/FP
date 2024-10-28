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

#################################
    #### Drehschieberpumpe

if DP == True:

    # t und D0 = usecols(2, 3)

    # Fehler für Evakuierung
    # 1000 - 10 mbar: 0.3%
    # 10 - 2*10^-3:    10%
    pE = 3.83E-3 # Enddruck
    VD = ufloat(34, 34*0.1) # Volumen Drehschieberpumpe

    t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Evak_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0] # Verschieb Zeit 
    D0_err = np.concatenate((D0[D0>=10]*0.003, D0[D0<10]*0.1)) # D0_err für Plots (nicht nützlich weil kaum erkennbar)


    ### Allgemeiner Plot 

    plt.errorbar(t, D0, yerr=D0_err, errorevery=(10), fmt=".", ecolor="pink", label="Messung")

    # Zeiten für LinFit: 0 - 150, 150 - 250, 250 - Ende 
    # 0 - 150
    params, cov = curve_fit(linfit, t[0:150], np.log((D0[0:150]-pE)/(D0[0]-pE)), sigma=D0_err[0:150], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x0150 = np.linspace(0, 170)
    plt.errorbar(x0150, pE + (D0[0] - pE)*np.exp(linfit(x0150, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
    print("DP Fit 1 S = ", S, " l/s")
    print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    print(" ")

    # 150 - 250
    params, cov = curve_fit(linfit, t[150:250], np.log((D0[150:250]-pE)/(D0[0]-pE)), sigma=D0_err[150:250], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x150250 = np.linspace(130, 300)
    plt.errorbar(x150250, pE + (D0[0] - pE)*np.exp(linfit(x150250, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
    print("DP Fit 2 S = ", S, " l/s")
    print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    print(" ")

    # 250 - Ende
    params, cov = curve_fit(linfit, t[250:], np.log((D0[250:]-pE)/(D0[0]-pE)), sigma=D0_err[250:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x250end = np.linspace(230, t[-1])
    plt.errorbar(x250end, pE + (D0[0] - pE)*np.exp(linfit(x250end, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
    print("DP Fit 3 S = ", S, " l/s")
    print("m = ", linfit_params[0], "   b = ", linfit_params[1])
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
    print("DP Leck 0.5 mbar S = ", S, " l/s")
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
    print("DP Leck 10 mbar S = ", S, " l/s")
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
    print("DP Leck 50 mbar S = ", S, " l/s")
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
    t1, D01 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck1_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t1 = t1 - t1[0]
    D01_err = D01*0.003
    D01 = unp.uarray(D01, D01_err)

    t2, D02 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck2_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t2 = t2 - t2[0]
    D02_err = D02*0.003
    D02 = unp.uarray(D02, D02_err)
    
    t3, D03 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck3_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t3 = t3 - t3[0]
    D03_err = D03*0.003
    D03 = unp.uarray(D03, D03_err)

    plt.errorbar(t1, unp.nominal_values(D01), yerr=D01_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung 1")
    plt.errorbar(t2, unp.nominal_values(D02), yerr=D02_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung 2")
    plt.errorbar(t3, unp.nominal_values(D03), yerr=D03_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung 3")
    plt.xlabel("t in s")
    plt.ylabel("p in mbar")
    plt.legend()
    plt.savefig("DP_Leck_100mbar_alle.pdf")
    plt.figure()

    # Mittelwertbildung 100 mbar
    N = min(len(t1), len(t2), len(t3))
    D01 = D01[:N]
    D02 = D02[:N]
    D03 = D03[:N]
    mean = (unp.nominal_values(D01) + unp.nominal_values(D02) + unp.nominal_values(D03))/3
    std_all = np.array((unp.std_devs(D01), unp.std_devs(D02), unp.std_devs(D03)))
    std_mean = np.sqrt(np.sum(std_all**2, axis=0)/N)
    t = t1[:N]
    D0 = unp.uarray(mean, std_mean)

        # plot für 100 mbar mittelwerte
    plt.errorbar(t, unp.nominal_values(D0), yerr=unp.std_devs(D0), fmt=".", ecolor="pink", errorevery=(5), label="Mittelwert")

        # fit für 100 mbar mittelwerte
    params, cov = curve_fit(linfit, t[10:150], unp.nominal_values(D0[10:150]), sigma=unp.std_devs(D0[10:150]), absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VD/pG * linfit_params[0]
    print("DP Leck Mittel, 100 mbar S = ", S, " l/s")
    print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    print(" ")

    plt.xlabel("Zeit t in s")
    plt.ylabel("Druck p in mbar")
    plt.title(f"p0 = {D0[0]} mbar, pE = {D0[-1]} mbar")
    plt.legend()
    plt.savefig("DP_Leck1_100mbar_mittelwert.pdf")
    plt.figure()



'''
DP Fit 1 S =  1.39+/-0.14  l/s
m =  -0.0408+/-0.0007    b =  0.01+/-0.07
 
DP Fit 2 S =  0.49+/-0.05  l/s
m =  -0.01441+/-0.00034    b =  -4.16+/-0.08
 
DP Fit 3 S =  0.41+/-0.04  l/s
m =  -0.011943+/-0.000008    b =  -4.116+/-0.004
 
DP Fit 4 S =  0.203+/-0.020  l/s
m =  -0.005971+/-0.000004    b =  -6.174+/-0.006
 
 
DP Leck 0.5 mbar S =  1.02+/-0.15  l/s
m =  0.0150+/-0.0004    b =  1.63+/-0.04
 
DP Leck 10 mbar S =  1.60+/-0.17  l/s
m =  0.47184+/-0.00022    b =  14.866+/-0.016
 
DP Leck 50 mbar S =  1.64+/-0.17  l/s
m =  2.4190+/-0.0010    b =  59.78+/-0.07
 
DP Leck Mittel, 100 mbar S =  1.54+/-0.16  l/s
m =  4.52830+/-0.00032    b =  112.640+/-0.018
'''
###################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
# 100 - 10^-8 mbar: 30%

pE = 4.0E-6 # Enddruck
VT = ufloat(33, 33*0.1) # Volumen Drehschieberpumpe

t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak1_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0] # Verschieb Zeit 
D2_err = D2*0.3


### Allgemeiner Plot 

plt.errorbar(t, D2, yerr=D2_err, errorevery=(5), fmt=".", ecolor="pink", label="Messung")

# Zeiten für LinFit: 0 - 15, 15 - 50, 50 - Ende 
# 0 - 15


params, cov = curve_fit(linfit, t[0:15], np.log((D2[:15]-pE)/(D2[0]-pE)), sigma=D2_err[0:15], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x015 = np.linspace(0, 20)
plt.errorbar(x015, pE + (D2[0] - pE)*np.exp(linfit(x015, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
print("TP Fit 1 S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 15 - 50
params, cov = curve_fit(linfit, t[15:50], np.log((D2[15:50]-pE)/(D2[0]-pE)), sigma=D2_err[15:50], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x1550 = np.linspace(10, 60)
plt.errorbar(x1550, pE + (D2[0] - pE)*np.exp(linfit(x1550, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
print("TP Fit 2 S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")

# 50 - Ende
params, cov = curve_fit(linfit, t[50:], np.log((D2[50:]-pE)/(D2[0]-pE)), sigma=D2_err[50:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
S = -linfit_params[0]*VT
x50 = np.linspace(40, t[-1])
plt.errorbar(x50, pE + (D2[0] - pE)*np.exp(linfit(x50, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
print("TP Fit 3 S = ", S, " l/s")
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
print("TP Leck 1E-4 mbar S = ", S, " l/s")
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
print("TP Leck 2E-4 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Leck_2e4.pdf")
plt.figure()

# # 5e-5 mbar
pG = unp.uarray(5E-5, (5E-5)*0.3)
t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_5e5_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D2_err = D2*0.3
plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VT/pG * linfit_params[0]
print("TP Leck 5E-5 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Leck_5e5.pdf")
plt.figure()

# # 7e-5 mbar
pG = unp.uarray(7E-5, (7E-5)*0.3)
t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_7e5_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
t = t - t[0]
D2_err = D2*0.3
plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="pink", errorevery=(5), label="Messung")

params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
x_fit = np.linspace(0, t[-1]+10)
plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
S = VT/pG * linfit_params[0]
print("TP Leck 7E-5 mbar S = ", S, " l/s")
print("m = ", linfit_params[0], "   b = ", linfit_params[1])
print(" ")
plt.xlabel("Zeit t in s")
plt.ylabel("Druck p in mbar")
plt.title(f"p0 = {D2[0]} mbar, pE = {D2[-1]} mbar")
plt.legend()
plt.savefig("TP_Leck_7e5.pdf")
plt.figure()

### Leitwerte
# t, D1, D2 = usecols(2, 4, 5)

# Fehler für Evakuierung:
# 100 - 10^-8 mbar: 30%
