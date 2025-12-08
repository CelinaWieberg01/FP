import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

deg0, t0, c0 = np.genfromtxt("data/data0.txt", unpack=True)
C0 = c0/t0
C0_err = np.sqrt(c0)/t0

degAu, tAu, cAu = np.genfromtxt("data/data_Au.txt", unpack=True)
CAu = cAu/tAu
CAu_err = np.sqrt(cAu)/tAu

degBi, tBi, cBi = np.genfromtxt("data/data_Bi.txt", unpack=True)
CBi = cBi/tBi
CBi_err = np.sqrt(cBi)/tBi

plt.errorbar(deg0, C0,  yerr=C0_err, fmt="o", c="blue", label='Vakuum', markersize=4)
plt.plot(deg0, C0, c="blue", alpha=0.15)

plt.xlabel('Winkel theta in °')
plt.ylabel('Zählrate C in 1/s')
plt.title('Zählrate in Abhängigkeit des Winkels, Vakuum')
plt.grid()
plt.legend()
plt.savefig('plots/plot0.pdf')


plt.figure()

plt.errorbar(degAu, CAu,  yerr=CAu_err, fmt="o", c="red", label='Goldfolie', markersize=4)
plt.plot(degAu, CAu, c="red", alpha=0.15)

plt.xlabel('Winkel theta in °')
plt.ylabel('Zählrate C in 1/s')
plt.title('Zählrate in Abhängigkeit des Winkels, Goldfolie')
plt.grid()
plt.legend()
plt.savefig('plots/plot_Au.pdf')

plt.figure()

plt.errorbar(degBi, CBi,  yerr=CBi_err, fmt="o", c="green", label='Bismutfolie', markersize=4)
plt.plot(degBi, CBi, c="green", alpha=0.15)
plt.xlabel('Winkel theta in °')
plt.ylabel('Zählrate C in 1/s')
plt.title('Zählrate in Abhängigkeit des Winkels, Bismutfolie')
plt.grid()
plt.legend()
plt.savefig('plots/plot_Bi.pdf')