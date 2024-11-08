import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
from scipy.signal import find_peaks

def impulserror(impulse):
    error = np.sqrt(impulse)
    array = unp.uarray(impulse, error)
    return array

# Europium
# Ganze Messung
Eu = np.genfromtxt("Europium.Spe", unpack=True)
Eu_err = impulserror(Eu)
channels = np.arange(len(Eu))

plt.plot(channels, Eu, label="Messwerte 152-Eu")
plt.plot([1187, 2379, 3343, 3985, 4307, 7549], [Eu[1187], Eu[2379], Eu[3343], Eu[3985], Eu[4307], Eu[7549]], "x", color="red", label="Ausgewählte Peaks")
FPE_index = 1187
plt.title("Messung für 152-Europium")
plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.yscale("log")
plt.legend()
plt.savefig("plots/Europium.pdf")
plt.figure()

# Einzelner Impuls
peak_channels = channels[FPE_index-35: FPE_index+35]
peak_impulses = Eu[FPE_index-35:FPE_index+35]
plt.plot(peak_channels, peak_impulses, "x", markersize=4, label="Peak")

untergrund_indices = np.where(peak_impulses < 40)[0]
untergrundkanal = peak_channels[untergrund_indices]
untergrundimpuls = peak_impulses[untergrund_indices]
plt.plot(untergrundkanal, untergrundimpuls, "x", markersize=4, color="red", label="Untergrund")


untergrund = impulserror(untergrundimpuls)
untergrund_mittel = np.sum(untergrund)/len(untergrund)
print(untergrund_mittel)
# untergrund_mittel = 21.5 +- 0.8

plt.plot()

signal_indices = np.where(peak_impulses > 40)[0]
signalkanal = peak_channels[signal_indices]
signalimpuls = peak_impulses[signal_indices]
plt.plot(signalkanal, signalimpuls - unp.nominal_values(untergrund_mittel), "x", markersize=4, color="green", label="Peak, bereinigt")

plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.legend()
plt.title("Peak bei 121.78 keV")
plt.savefig("plots/PeakUntergrund.pdf")
plt.figure()

# Linieninhalte
peaks = [1187, 2379, 3343, 3985, 4307, 7549]
linieninhalte_Eu = unp.uarray(np.zeros(6), np.zeros(6))
for i, peak in enumerate(peaks):
    linieninhalte_Eu[i] = sum(Eu_err[peak-20:peak+20])

peaks_Energie = [121.78, 244.70, 344.30, 411.12, 443.96, 778.90]

plt.plot(peaks, peaks_Energie, "x", c="red", label="Peak-Energien")
params = np.polyfit(peaks, peaks_Energie, 1)
xx = np.linspace(0, 8000, 1000)
plt.plot(xx, params[0]*xx + params[1], c="green", label=f"Fit mit m = {params[0]}, b = {params[1]}")
plt.xlabel("Kanal")
plt.ylabel("Energie")
plt.legend()
plt.savefig("plots/Kalibrierung.pdf")
plt.figure()

"""
Energie (keV)   Wahrschl (%) Kanalnummer Linieninhalt
121.78   28.6  1187  10375 +- 101
244.70   7.6   2379  1905 +- 43
344.30   26.5  3343  4081 +- 63
411.12   2.2   3985  406 +- 20
443.96   3.1   4307  488 +- 22
778.90   12.9  7549  841 +- 29
"""
W_Eu = [28.6, 7.6, 26.5, 2.2, 3.1, 12.9]
# Aktivität
t_h_152Eu = 13 + 196/365
t_datum = 24.11
A_heute = ufloat(4130, 60) * np.exp(-np.log(2)*t_datum/t_h_152Eu)
"""
Aktivität = 1202 +- 17 Bq
"""
abstand = ufloat(70.2, 0.05)
Omega = 2*np.pi*(1 - (15 + abstand)/(unp.sqrt((15+abstand)**2 + 22.5**2)))
"""
Omega = 0.20826 +- 0.00023
"""
t_mess = 3553
Q = np.array(linieninhalte_Eu)*4*np.pi/(np.array(W_Eu) * A_heute * t_mess * Omega)
"""
Q = [0.005126487218972472+/-9.007010076499762e-05
 0.003542249881436617+/-9.617962610705257e-05
 0.0021762988247944447+/-4.6541345494713615e-05
 0.0026079614016131765+/-0.00013489383611213702
 0.0022246204995804123+/-0.0001057919690081008
 0.0009213064087647879+/-3.448892266253558e-05]
"""
def effizienz(E, alpha, beta):
    return alpha*E**(beta)

params, cov = curve_fit(effizienz, peaks_Energie, unp.nominal_values(Q), sigma=unp.std_devs(Q), absolute_sigma=True)
print(f"alpha = {params[0]} +- {np.sqrt(np.diag(cov))[0]}, beta = {params[1]} +- {np.sqrt(np.diag(cov))[1]}")
xx = np.linspace(100, 800, 1000)
plt.plot(xx, effizienz(xx, params[0], params[1]), c="red", label="fit")
"""
alpha = 0.32903965480974107 +- 0.03072936509917969, beta = -0.8572976119270979 +- 0.0169677843756542
"""
plt.errorbar(peaks_Energie, unp.nominal_values(Q), yerr=unp.std_devs(Q), fmt="x", label="Effizienz")
plt.xlabel("Energie E in keV")
plt.ylabel("Effizienz Q")
plt.legend()
plt.savefig("plots/Effizienz.pdf")