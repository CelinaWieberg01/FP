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

def linfit(x, m, b):
    return m*x + b

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
#print(untergrund_mittel)
# untergrund_mittel = 21.5 +- 0.8

plt.plot()

signal_indices = np.where(peak_impulses > 40)[0]
signalkanal = peak_channels[signal_indices]
signalimpuls = peak_impulses[signal_indices]
plt.plot(signalkanal, signalimpuls - unp.nominal_values(untergrund_mittel), "x", markersize=4, color="green", label="Peak, bereinigt")

plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.legend()
plt.title("Eu-Peak bei 121.78 keV")
plt.savefig("plots/EuPeakUntergrund.pdf")
plt.figure()

# Linieninhalte
peaks = [1187, 2379, 3343, 3985, 4307, 7549]
linieninhalte_Eu = unp.uarray(np.zeros(6), np.zeros(6))
for i, peak in enumerate(peaks):
    linieninhalte_Eu[i] = sum(Eu_err[peak-20:peak+20])

peaks_Energie = [121.78, 244.70, 344.30, 411.12, 443.96, 778.90]

plt.plot(peaks, peaks_Energie, "x", c="red", label="Peak-Energien")
params_linfit, cov_linfit  = curve_fit(linfit, peaks, peaks_Energie, )
xx = np.linspace(0, 8000, 1000)
plt.plot(xx, params_linfit[0]*xx + params_linfit[1], c="green")
#print(f"Fit mit m = {params_linfit[0]} +- {np.sqrt(np.diag(cov))[0]}, b = {params_linfit[1]} +- {np.sqrt(np.diag(cov))[1]}")
"""
Fit mit m = 0.10330462983489719 +- 4.407553331752846e-05, b = -0.9033881314038474 +- 0.18842996495320607
"""
plt.xlabel("Kanal")
plt.ylabel("Energie")
plt.title("Kalibrierungskurve")
plt.legend()
plt.savefig("plots/EuKalibrierung.pdf")
plt.figure()

def kanalzuenergie(kanal):
    return ufloat(params_linfit[0], np.sqrt(np.diag(cov_linfit))[0]) * kanal + ufloat(params_linfit[1], np.sqrt(np.diag(cov_linfit))[1])

"""
Energie (keV)   Wahrschl (%) Kanalnummer Linieninhalt
121.78   28.6  1188  10375 +- 101
244.70   7.6   2380  1905 +- 43
344.30   26.5  3344  4081 +- 63
411.12   2.2   3986  406 +- 20
443.96   3.1   4308  488 +- 22
778.90   12.9  7550  841 +- 29
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

params_eff, cov_eff = curve_fit(effizienz, peaks_Energie, unp.nominal_values(Q), sigma=unp.std_devs(Q), absolute_sigma=True)
#print(f"alpha = {params[0]} +- {np.sqrt(np.diag(cov))[0]}, beta = {params[1]} +- {np.sqrt(np.diag(cov))[1]}")
xx = np.linspace(100, 800, 1000)
plt.plot(xx, effizienz(xx, params_eff[0], params_eff[1]), c="red", label="fit")
"""
alpha = 0.32903965480974107 +- 0.03072936509917969, beta = -0.8572976119270979 +- 0.0169677843756542
"""
plt.errorbar(peaks_Energie, unp.nominal_values(Q), yerr=unp.std_devs(Q), fmt="x", label="Effizienz")
plt.xlabel("Energie E in keV")
plt.ylabel("Effizienz Q")
plt.title("Effizienzkurve")
plt.legend()
plt.savefig("plots/EuEffizienz.pdf")
plt.figure()

#### 
# 137-Cs

Cs = np.genfromtxt("Caesium.Spe", unpack=True)
Cs_err = impulserror(Cs)
channels = np.arange(len(Cs))
Cs_err = impulserror(Cs)

plt.plot(channels, Cs, label="Messwerte 137-Cs")
plt.plot([4430], Cs[4430], "x", color="red", label="Comptonkante")
plt.plot([1856], Cs[1856], "x", color="green", label="Rückstreupeak")
plt.plot([6414], Cs[6414], "x", color="orange", label="Photopeak")
plt.title("Messung für 137-Cäsium")
plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.yscale("log")
plt.legend()
plt.savefig("plots/Caesium.pdf")
plt.figure()

# peak bei channels[6414] also in Kanal 6415
channels_peak = channels[6414-45:6414+45]
Cs_peak = Cs[6414-45:6414+45]
Cs_peak_err = impulserror(Cs_peak)
Cs_peak_err[Cs_peak_err == ufloat(0, 0)] = ufloat(0,1)
plt.plot(channels_peak, Cs_peak, "x", label="Messwerte")

Cs_peak_inhalt = np.sum(Cs_peak_err[25:-25])
#print("Cs_peak_inhalt = ", Cs_peak_inhalt)
"""
Cs_peak_inhalt =  (1.176+/-0.011)e+04
"""

def gauss(x, a, my, sigma):
    return a/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x - my)**2 / (2*sigma**2))

params, cov = curve_fit(gauss, channels_peak, Cs_peak, p0=[13000, 6417, 10], sigma=unp.std_devs(Cs_peak_err), absolute_sigma=True)
#print(f"a = {params[0]} +- {np.sqrt(np.diag(cov))[0]}, my = {params[1]} +- {np.sqrt(np.diag(cov))[1]}, sigma = {params[2]} +- {np.sqrt(np.diag(cov))[2]}")
"""
a = 12233.18358505073 +- 110.63744277561693, my = 6415.985582417949 +- 0.08947257593154681, sigma = 9.764919544935347 +- 0.06341548350152967
"""
xx = np.linspace(6360, 6460, 1000)
fit_y = gauss(xx, *params)
plt.plot(xx, fit_y, color="red", label="Fit")

Cs_max = max(fit_y)
Cs2 = Cs_max/2
#print("Cs2 = ", Cs2)
Cs10 = Cs_max/10
#print("Cs10 = ", Cs10)

channels_Cs2 = xx[fit_y > Cs2]
plt.hlines(Cs2, channels_Cs2[0], channels_Cs2[-1], color="orange", label="Halbwertsbreite")
channels_Cs10 = xx[fit_y > Cs10]
plt.hlines(Cs10, channels_Cs10[0], channels_Cs10[-1], color="green", label="Zehntelwertsbreite")

Cs2_breite = channels_Cs2[-1] - channels_Cs2[0]
#print("Cs2_breite = ", Cs2_breite)
Cs2_energie = kanalzuenergie(Cs2_breite)
#print("Cs2_energie = ", unp.nominal_values(Cs2_energie), " +- ", unp.std_devs(Cs2_energie))

Cs10_breite = channels_Cs10[-1] - channels_Cs10[0]
#print("Cs10_breite = ", Cs10_breite)
Cs10_energie = kanalzuenergie(Cs10_breite)
#print("Cs10_energie = ", unp.nominal_values(Cs10_energie), " +- ", unp.std_devs(Cs10_energie))
#print("Checking for plausability. E10 / E2 = ", Cs10_energie/Cs2_energie)
"""
Cs2 =  249.89001045659492
Cs10 =  49.978002091318984
Cs2_breite =  22.922922922922226
Cs2_energie =  1.4646559358825124  +-  0.18843267359799684
Cs10_breite =  41.841841841841415
Cs10_energie =  3.4190678516778927  +-  0.1884389895202624
Checking for plausability. E10 / E2 =  2.33+/-0.33
"""

plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.title("Cs-Peak bei 662 keV")
plt.legend()
plt.savefig("plots/CsPeak.pdf")
plt.figure()

# Compton-Kante und Rückstreupeak
FEP_Cs_err = impulserror(6414)
E_FEP_Cs = kanalzuenergie(FEP_Cs_err)
#print("E_FEP_Cs = ", E_FEP_Cs)

Compton_Cs_err = impulserror(4430)
E_Compton_Cs = kanalzuenergie(Compton_Cs_err)
#print("E_Compton_Cs = ", E_Compton_Cs)

Back_Cs_err = impulserror(1856)
E_Back_Cs = kanalzuenergie(Back_Cs_err)
#print("E_Back_Cs = ", E_Back_Cs)
"""
E_FEP_Cs =  662+/-8
E_Compton_Cs =  457+/-7
E_Back_Cs =  191+/-4
"""

m_e = constants.electron_mass
c = constants.speed_of_light
E_m = m_e*c**2
epsilon = FEP_Cs_err / E_m

Theorie_Compton_Kante = E_FEP_Cs * 2 * epsilon / (1 + 2 * epsilon)
#print("Theorie_Compton_Kante = ", Theorie_Compton_Kante)

Theorie_Back = E_FEP_Cs / (1 + 2*epsilon)
#print("Theorie_Back = ", Theorie_Back)

channels_diffwirk = channels[2000:4000]
Cs_diffwirk = Cs_err[2000:4000]
energie_channels = kanalzuenergie(channels_diffwirk)

def querschnitt(E, k):
    dsigmadE = k*(2 + (E / (FEP_Cs_err - E))**2 * (1/epsilon**2 + (FEP_Cs_err - E)/(FEP_Cs_err) - (2/epsilon)*((FEP_Cs_err - E)/E)))
    return dsigmadE

#params, cov = curve_fit(querschnitt, unp.nominal_values(energie_channels), unp.nominal_values(Cs_diffwirk))
# Fit nicht möglich boah ich bring mich um :(

# Absorptionswahrscheinlichkeiten

def abswsl(my, l):
    return 1 - np.exp(-my*l)

my_P = 0.007
my_C = 0.37
d = 3.9
#print("Wsl_Photo = ", abswsl(my_P, d))
#print("Wsl_Compton = ", abswsl(my_C, d))

########

# Barium-133

Ba = np.genfromtxt("Barium.Spe", unpack=True)
Ba_err = impulserror(Ba)
channels = np.arange(len(Ba))
t = 3461
Ba_peaks = [793, 2684, 2944, 3456, 3727]
Ba793 = kanalzuenergie(793)
Ba2684 = kanalzuenergie(2684)
Ba2944 = kanalzuenergie(2944)
Ba3456 = kanalzuenergie(3456)
Ba3727 = kanalzuenergie(3727)

Ba_energien = (Ba793, Ba2684, Ba2944, Ba3456, Ba3727)

Linieninhalte = np.empty(5, dtype=type(Ba_err))

for i, peak in enumerate(Ba_peaks):
    Linieninhalt = np.sum(Ba_err[peak-19:peak+21])
    Linieninhalte[i] = Linieninhalt
#print(Linieninhalte)

Wahrscheinlichkeiten = unp.uarray((33.31, 7.13, 18.31, 62.05, 8.94), (0.30, 0.06, 0.11, 0.19, 0.06))
#print(Wahrscheinlichkeiten)

Effizienz = effizienz(Ba_energien, *params_eff)
#print(Effizienz)

def Aktivitaet(Q, W, Z):
    return 4*np.pi*Z/(Q*W*t*Omega)

A_Ba = np.mean(Aktivitaet(Effizienz, Wahrscheinlichkeiten, Linieninhalte))
#print(A_Ba)

plt.plot(channels, Ba, label="Barium-133")
plt.plot(Ba_peaks, Ba[Ba_peaks], "x", color="red", label="Peaks")

plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.title("Messung für Barium-133")
plt.yscale("log")
plt.legend()
plt.savefig("plots/Barium.pdf")
plt.figure()

######################
# Uranophan

H = np.genfromtxt("Hauistdumm.Spe", unpack=True)
H_err = impulserror(H)
channels=np.arange(len(H))

plt.plot(channels, H, label="Mineral")

H_peaks = np.array((755, 903, 1808, 2351, 2867, 3417, 5913, 6455, 7442, 7626, 7818))
plt.plot(H_peaks, H[H_peaks], "x", color="red", label="Peaks")

H_peaks_energien = kanalzuenergie(H_peaks)
#print(H_peaks_energien)

Linieninhalte = np.empty(len(H_peaks), dtype=type(Ba_err))

for i, peak in enumerate(H_peaks):
    Linieninhalt = np.sum(H_err[peak-19:peak+21])
    Linieninhalte[i] = Linieninhalt
#print(Linieninhalte)

Wahrscheinlichkeiten = unp.uarray(np.array((3.75, 2.18, 3.555, 7.268, 18.414, 35.6, 45.49, 1.530, 4.892, 1.064, 1.262)), np.array((0.08, 0.19, 0.019, 0.022, 0.036, 0.7, 0.19, 0.007, 0.016, 0.013, 0.006)))
#print(Wahrscheinlichkeiten)

t = 3582

Effizienz = effizienz(H_peaks_energien, *params_eff)
#print(Effizienz)

A_H = (Aktivitaet(Effizienz, Wahrscheinlichkeiten, Linieninhalte))

#print(A_H)



plt.xlabel("Kanal")
plt.ylabel("Impulse")
plt.title("Messung für unbekanntes Erz")
plt.yscale("log")
plt.legend()
plt.savefig("plots/H.pdf")