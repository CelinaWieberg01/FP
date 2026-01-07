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
print("gauss: f(x) = a * e^((x-x0)^2 / (2*sigma^2)")
print(f"a = {a:.10} Einheit: Counts")
print(f"x0 = {x0:.10} Einheit: ns")
print(f"sigma = {sigma:.10} Einheit: ns")

fwhm = 2*np.sqrt(2*np.log(2)) * (sigma)
print(f"fwhm breite = {fwhm:.10} Einheit: ns")
print(" ")
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

m = ufloat(p_cal[0], np.sqrt(np.diag(cov_cal))[0])
b = ufloat(p_cal[1], np.sqrt(np.diag(cov_cal))[1])
print("linfit für kalibrierung: t(K) = m*K + b")
print("m = ", m, " Einheit: µs/Kanalnummer")
print("b = ", b, " Einheit: µs")
print(" ")

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




counts = np.genfromtxt("Daten.Spe", unpack=True)
channels = np.arange(1, len(counts)+1, 1)
zeiten = linfit(channels, *p_cal)

def zerfall(x, a, b, c, U):
    return a*np.exp(-b*(x-c)) + U

params_e, cov_e = curve_fit(zerfall, zeiten[4:400], counts[4:400], p0=(40, 0.5, 0, 0))

a = ufloat(params_e[0], np.sqrt(np.diag(cov_e))[0])
b = ufloat(params_e[1], np.sqrt(np.diag(cov_e))[1])
c = ufloat(params_e[2], np.sqrt(np.diag(cov_e))[2])
U = ufloat(params_e[3], np.sqrt(np.diag(cov_e))[3])

print("zerfallskurve = a*e^(-lambda*t - c)")
print(f"a =  {a:.10}, Einheit: Counts")
print(f"lambda = {b:10} Einheit: 1/µs")
print(f"c = {c:.10} Einheit: µs")
print(f"U = {U:.10} Einheit: Counts")
print(" ")
print(f"Mittlere Lebensdauer = 1/lambda = {1/b :.10} Einheit: µs")


xx = np.linspace(0, max(zeiten))
yy = zerfall(xx, *params_e)
plt.plot(xx, yy, color="black", label="Fit", zorder=10)

plt.errorbar(
    zeiten[:4],
    counts[:4],
    yerr=np.sqrt(counts[:4]),
    fmt=".",
    markersize=5,
    color="red",
    ecolor=(1, 0, 0, 0.3),  # transparent red error bars
    label="Unberücksichtigt"
)

plt.errorbar(
    zeiten[4:400],
    counts[4:400],
    yerr=np.sqrt(counts[4:400]),
    fmt=".",
    markersize=5,
    ecolor=(0, 0, 1, 0.3),  # transparent black error bars
    label="Berücksichtigt"
)

plt.errorbar(
    zeiten[400:],
    counts[400:],
    yerr=np.sqrt(counts[400:]),
    fmt=".",
    markersize=5,
    color="red",
    ecolor=(1, 0, 0, 0.3)
)


plt.xlabel(r"Zeit in \si{\micro\second}")
plt.ylabel("Counts")
plt.grid("on")
plt.title("Lebensdauer einzelner Myonen")
plt.legend()

plt.savefig("plots/zerfall.pdf")

N_Myon = 1496787/175057
T_s = 1e-5
P1 = T_s * N_Myon * np.exp(T_s * N_Myon)
U = 1496787*P1/512
print(f"U = {U} Einheit: Counts/Kanal")

"""
gauss: f(x) = a * e^((x-x0)^2 / (2*sigma^2)
a = 283.6143354+/-7.8672758 Einheit: Counts
x0 = 0.6513704308+/-0.1750581592 Einheit: ns
sigma = 5.066269628+/-0.196916370 Einheit: ns
fwhm breite = 11.93015327+/-0.46370262 Einheit: ns
 
linfit für kalibrierung: t(K) = m*K + b
m =  0.021700+/-0.000015  Einheit: µs/Kanalnummer
b =  0.1577+/-0.0035  Einheit: µs
 
zerfallskurve = a*e^(-lambda*t - c)
a =  99.9+/-160826667.0, Einheit: Counts
lambda =      0.543+/-     0.020 Einheit: 1/µs
c = -1.981+/-2963191.422 Einheit: µs
 
Mittlere Lebensdauer = 1/lambda = 1.840427112+/-0.068478126 Einheit: µs
"""