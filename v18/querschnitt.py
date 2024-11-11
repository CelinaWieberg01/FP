import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl


def impulserror(impulse):
    error = np.sqrt(impulse)
    array = unp.uarray(impulse, error)
    return array

def kanalzuenergie(kanal):
    return ufloat(0.10330462983489719, 4.407553331752846e-05) * kanal + ufloat(-0.9033881314038474, 0.18842996495320607)


Cs = np.genfromtxt("Caesium.Spe", unpack=True)
Cs_err = impulserror(Cs)
channels = np.arange(len(Cs))
E_channels = kanalzuenergie(channels)

plt.plot(unp.nominal_values(E_channels), Cs, label="Messungen")


# FEP
FEP_Cs_err = impulserror(6414)
E_FEP_Cs_err = unp.nominal_values(kanalzuenergie(FEP_Cs_err))
print("E_FEP_Cs = ", E_FEP_Cs_err)

# Comptonkante
Compton_Cs_err = impulserror(4430)
E_Compton_Cs_err = kanalzuenergie(Compton_Cs_err)
print("E_Compton_Cs = ", E_Compton_Cs_err)

# Rückstreupeak
Back_Cs_err = impulserror(1856)
E_Back_Cs_err = kanalzuenergie(Back_Cs_err)
print("E_Back_Cs = ", E_Back_Cs_err)

# Werte für theoretische Messung
#e_charge = constants.elementary_charge
m_e = 511
c = constants.speed_of_light
E_m = m_e
epsilon = unp.nominal_values(E_FEP_Cs_err / E_m)

# Theoretischer Wert für Compton-Kante
Theorie_Compton_Kante = E_FEP_Cs_err * 2 * epsilon / (1 + 2 * epsilon)
print("Theorie_Compton_Kante = ", Theorie_Compton_Kante)

# Theoretischer Wert für Rückstreupeak
Theorie_Back = E_FEP_Cs_err / (1 + 2*epsilon)
print("Theorie_Back = ", Theorie_Back)

# Kanäle für Fit
channels_diffwirk = channels[2300:4000]
energie_channels = kanalzuenergie(channels_diffwirk)

Cs_diffwirk_err = Cs_err[2300:4000]

def querschnitt(E, k):
    dsigmadE = k*(2 + (E / (E_FEP_Cs_err - E))**2 * (1/epsilon**2 + (E_FEP_Cs_err - E)/(E_FEP_Cs_err) - (2/epsilon)*((E_FEP_Cs_err - E)/E)))
    return dsigmadE

params, cov = curve_fit(querschnitt, unp.nominal_values(energie_channels), unp.nominal_values(Cs_diffwirk_err))

xx = np.linspace(0, 4600, 1000)
xx = unp.nominal_values(kanalzuenergie(xx))
plt.plot(xx, querschnitt(xx, *params), color="red", label="Fit")
plt.xlabel(r"Energie in $\si{\kilo\electronvolt}$")
plt.ylabel("Impulse")
plt.yscale("log")
plt.legend()
plt.savefig("plots/querschnitt.pdf")
print(np.sum(querschnitt(xx, *params)))
print(params[0], " +- ", np.sqrt(cov))