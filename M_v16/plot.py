import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

background = False
if background == True:
    bgAu = 0.31
    bgBi = 0.54
else:
    bgAu = bgBi = 0

deg0, t0, c0 = np.genfromtxt("data/data0.txt", unpack=True)
C0 = c0/t0
C0_err = np.sqrt(c0)/t0

def normal(x, a, b):
    return a*np.exp(-b*(x)**2)

mask0 = [6,7,8,10,12,13,14]
deg0_fit = deg0[mask0]
C0_fit   = C0[mask0]

params0, cov0 = curve_fit(normal, deg0_fit, C0_fit)
theta_fit = np.linspace(-4, 9, 1000)
C_fit = normal(theta_fit, *params0)

print(f"FIT CHARAKTERISTIK \na = {params0[0]} +- {np.sqrt(np.diag(cov0))[0]} \nb = {params0[1]} +- {np.sqrt(np.diag(cov0))[1]}\n")

def bg(theta):
    return normal(theta, *params0)

degAu, tAu, cAu = np.genfromtxt("data/data_Au.txt", unpack=True)
CAu = cAu/tAu
CAu_err = np.sqrt(cAu)/tAu
CAu_werr = unp.uarray(CAu-bgAu*bg(degAu), CAu_err)

degBi, tBi, cBi = np.genfromtxt("data/data_Bi.txt", unpack=True)
CBi = cBi/tBi
CBi_err = np.sqrt(cBi)/tBi
CBi_werr = unp.uarray(CBi-bgBi*bg(degBi), CBi_err)

plt.errorbar(deg0, C0, yerr=C0_err, fmt="o", c="blue", label='Vakuum', markersize=4)
plt.plot(deg0, C0, c="blue", alpha=0.15)
plt.scatter(deg0_fit, C0_fit, facecolors="none", edgecolors="red", s=80, label="für Fit berücksichtigter Wert")
plt.plot(theta_fit, C_fit, c="orange", label="Fit")

plt.xlabel(r'Winkel $\theta$ in \si{\degree}')
plt.ylabel(r'Zählrate $C$ in \si{\per\second}')
plt.title('Zählrate in Abhängigkeit des Winkels, Vakuum')
plt.grid()
plt.legend()
plt.savefig('plots/meas_Vac.pdf')


plt.figure()

plt.errorbar(degAu, CAu, yerr=CAu_err, fmt="o", c="red", label='Goldfolie', markersize=4)
plt.plot(degAu, CAu, c="red", alpha=0.15)

plt.xlabel('Winkel theta in °')
plt.ylabel('Zählrate C in 1/s')
plt.title('Zählrate in Abhängigkeit des Winkels, Goldfolie')
plt.grid()
plt.legend()
plt.savefig('plots/meas_Au.pdf')

plt.figure()

plt.errorbar(degBi, CBi,  yerr=CBi_err, fmt="o", c="green", label='Bismutfolie', markersize=4)
plt.plot(degBi, CBi, c="green", alpha=0.15)
plt.xlabel('Winkel theta in °')
plt.ylabel('Zählrate C in 1/s')
plt.title('Zählrate in Abhängigkeit des Winkels, Bismutfolie')
plt.grid()
plt.legend()
plt.savefig('plots/meas_Bi.pdf')

plt.figure()

C0 = 330000*np.exp(-11391/157680) * 20/(4*np.pi*56**2) # Bq
n_au = constants.Avogadro*19320/0.197
n_bi = constants.Avogadro*9790/0.209
d = 2e-6
Omega = 4*np.arctan(10/90)*np.arctan(2/90)

epsilon0 = 8.8541878188e-12
e0 = 1.602176634e-19

def dcs_meas(C, n):   # C = Zählrate im Winkel Theta, n = Target-Atomkonzentration, d = Foliendicke, Omega = Raumwinkel, C0 = Rate einfallende Teilchen
    dcs = C / (n*d*Omega*C0)
    return dcs


dcs_au = dcs_meas(CAu_werr, n_au)
dcs_bi = dcs_meas(CBi_werr, n_bi)


def rutherford(theta, z, Z, K):
    dcs_th = (z*Z*e0**2/(4*np.pi*epsilon0*4*K))**2 * 1/(np.sin(theta*np.pi/360))**4
    return dcs_th

def rtf_fit(theta, a, b):
    return a/np.sin(np.pi*(theta-b)/360)**4 

theta = np.linspace(5,22)
z = 2
Z_au = 79
Z_bi = 83
K = 8.78954e-13

print("verwendete Parameter: z = 2 \nZ_au = 79 \nZ_bi = 83 \nK = 8.78954e-13 J")

params_Au, cov_Au = curve_fit(rtf_fit, degAu, unp.nominal_values(dcs_au))
print(f"Fit Au: \na = {params_Au[0]} +- {np.sqrt(np.diag(cov_Au))[0]} \nb = {params_Au[1]} +- {np.sqrt(np.diag(cov_Au))[1]}\n")
params_Bi, cov_Bi = curve_fit(rtf_fit, degBi, unp.nominal_values(dcs_bi), p0=(1e-19, 1))
print(f"Fit Bi: \na = {params_Bi[0]} +- {np.sqrt(np.diag(cov_Bi))[0]} \nb = {params_Bi[1]} +- {np.sqrt(np.diag(cov_Bi))[1]}\n")

rtf_au = rutherford(theta, z, Z_au, K)
rtf_fit_au = rtf_fit(theta, *params_Au)

rtf_bi = rutherford(theta, z, Z_bi, K)
rtf_fit_bi = rtf_fit(theta, *params_Bi)


plt.plot(theta, rtf_fit_au, color="gold", alpha=0.5, label="Fit")
plt.plot(theta, rtf_au, color="black", alpha = 0.5, label="Theoriekurve")
plt.errorbar(degAu, unp.nominal_values(dcs_au), yerr=unp.std_devs(dcs_au), fmt="o", c="red", label='Goldfolie', markersize=4)
plt.plot(degAu, unp.nominal_values(dcs_au), color="red", alpha=0.15)

plt.xlabel(r'Winkel $\theta$ in \si{\degree}')
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ in \si{\meter\square\per\steradian}")
plt.title('Goldfolie')
plt.grid()
plt.legend()
plt.savefig('plots/dsdO_Au.pdf')

plt.figure()

plt.plot(theta, rtf_fit_bi, color="gold", alpha=0.5, label="Fit")
plt.plot(theta, rtf_bi, color="black", alpha = 0.5, label="Theoriekurve")
plt.errorbar(degBi, unp.nominal_values(dcs_bi), yerr=unp.std_devs(dcs_bi), fmt="o", c="green", label='Bismutfolie', markersize=4)
plt.plot(degBi, unp.nominal_values(dcs_bi), color="green", alpha=0.15)

plt.xlabel(r'Winkel $\theta$ in \si{\degree}')
plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}$ in \si{\meter\square\per\steradian}")
plt.title('Bismutfolie')
plt.grid()
plt.legend()
plt.savefig('plots/dsdO_Bi.pdf')