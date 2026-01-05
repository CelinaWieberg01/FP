import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
from uncertainties import nominal_value as nom

# B im Delay, messdauer = 20 s
B   = np.array((0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10, 12, 16)) # ns
C_B = np.array((262, 237, 267, 338, 266, 223, 210, 235, 206, 205, 184, 173, 203, 184, 153, 162, 97 ,  21,  77, 62, 16,  6)) # num

# A im Delay, messdauer = 20 s
A   = np.array((0,   0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 16.0)) # ns
C_A = np.array((290, 318, 324, 271, 246, 235, 232, 208, 189, 167, 133, 125, 105,  76,  42,  53,   46,   17,   18)) # num


def gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

x = np.concatenate((-A, B))
y = np.concatenate((C_A, C_B))

popt, pcov = curve_fit(gauss, x, y)

a = ufloat(popt[0], np.sqrt(np.diag(pcov))[0])
x0 = ufloat(popt[1], np.sqrt(np.diag(pcov))[1])
sigma = ufloat(popt[2], np.sqrt(np.diag(pcov))[2])
print("a, x0, sigma")
print(a, x0, sigma)

fwhm = 2*np.sqrt(2*np.log(2)) * (sigma)
print("fwhm breite = ", fwhm)
xx = np.linspace(-20, 20)


plt.scatter(B, C_B, color="red", label=r"Delay in B")
plt.scatter(-A, C_A, color="blue", label=r"Delay in A")

plt.plot(xx, gauss(xx, *popt), label=r"Gaußfit", color="green")

plt.vlines([nom(x0 - fwhm/2), nom(x0 + fwhm/2)], 0, nom(a), linestyles="dashed", color="lightgreen", label=r"Grenze FWHM")
plt.hlines(nom(a/2), nom(x0 - fwhm/2), nom(x0 + fwhm/2), linestyles="dashed", color="orange", label=r"Höhe FWHM")

plt.xlabel(r"Delay in \si{\nano\second}")
plt.xlim(-20, 20)
plt.ylabel(r"Counts")
plt.ylim(top=350)
plt.grid("on")
plt.title(r"Peakabdeckung bei verschiedenen Delays")

plt.legend()

plt.savefig("plots/delay.pdf")
plt.figure()


# doppelimpulsgenrator MCA kalibrierung

impulsabstand = 0.1*np.array((10, 15, 20,  25,  30,  35,  40,  45,  50,  55,  60,  65,  70,  75,  80,  85)) # 0.1 µs
kanalnummer   = np.array((39, 62, 85, 108, 131, 154, 177, 200, 223, 246, 269, 292, 315, 338, 362, 385)) # num

def linfit(x, m, b):
    return m*x+b 

p_cal, cov_cal = curve_fit(linfit, kanalnummer, impulsabstand)

print(p_cal)
print(np.sqrt(np.diag(cov_cal)))

xx = np.linspace(0, 400)
yy = linfit(xx, *p_cal)

plt.plot(xx, yy, color="pink", label=r"Linearer Fit")
plt.scatter(kanalnummer, impulsabstand, color="purple", label=r"Messung")

plt.xlabel(r"Kanalnummer")
plt.xlim(0, 400)
plt.ylim(0, 9)
plt.ylabel(r"Impulsabstand in \si{\micro\second}")
plt.grid("on")
plt.title(r"Kalibrierungskurve MCA")
plt.legend()

plt.savefig("plots/kali.pdf")
plt.figure()