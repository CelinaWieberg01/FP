import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
plt.switch_backend('TkAgg')

#Datensatz 1
t1, T1 , I1= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten1.txt", unpack=True)
t1_orig, T1_orig , I1_orig= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten1.txt", unpack=True)
I1_pA= I1_orig/10 #umrechnung auf pA
T1=T1+273.15 #umrechnung auf Kelvin

#Linearer Fit für den ersten Untergrund
fit_T1 = np.concatenate((T1[12:17], T1[-16:-9]))  # Anfangs- und Endwerte für den Fit
fit_I1 = np.concatenate((I1[12:17], I1[-16:-9]))
fit_params1, fit_cov1 = np.polyfit(fit_T1, fit_I1, 1, cov=True)
alpha1 = ufloat(fit_params1[0], np.sqrt(fit_cov1[0][0]))
beta1 = ufloat(fit_params1[1], np.sqrt(fit_cov1[1][1]))
fit_line1_x = np.linspace(min(T1), max(T1), 100)
fit_line1_y = np.polyval(fit_params1, fit_line1_x)

# Linearer Fit für die erste Heizrate
heizrate1_parasm= np.polyfit(t1, T1, 1)  
heizrate1 = heizrate1_parasm[0]

#Original Plot für die erste Messreihe
plt.scatter(T1, I1, color='purple', marker='x', label='Erste Heizrate')
plt.plot(fit_line1_x, fit_line1_y, color='green', linestyle='-', label='Untergrund Fit') 
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Der Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 300)
plt.ylim(0, 5)
#plt.show()
plt.clf()

#Angabe Heizrate 1
print(f'Heizrate für die erste Messung: {alpha1.nominal_value:.2e} K/min')

#Datenbereinigung Reihe 1
untergrund1 = alpha1.nominal_value * T1 + beta1.nominal_value
bereinigt1 = I1-untergrund1
print(f'Erste Messung: alpha1 = {alpha1:.2e} pA/K, beta1 = {beta1:.2e} A')

#Plot der bereinigten Daten für Reihe 1
plt.scatter(T1, bereinigt1, color='purple', marker='x', label='Bereinigter erster Datensatz durch abziehen des Untergrunds')
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Bereinigter Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 340)
plt.ylim(0, max(bereinigt1))
#plt.show()
plt.clf()

#Aktivierungsenergie über den Polarisationsansatz
k_B = constants.Boltzmann

#Berechnung des Logarithmus von I und der Inversen Temperatur für Messreihe 1
I1_A= I1_orig * 1e-11#auf Ampere
ln_I1 = np.log(I1)
inv_T1= 1 /k_B * T1 

#Ausgabe
print("ln(I1):",ln_I1)
print("inv_T1(I1):", inv_T1)

# Plot für ln(I1) gegen 1/(kB*T1) ohne Ausgleichsgerade
plt.scatter(inv_T1, ln_I1, color='purple', marker='x', label='ln(I_1) gegen 1 / (k_B T)')
plt.xlabel(r'$1 / (k_B T) 1/J$')
plt.ylabel(r'$ln(I)$ (Einheitenlos)')
plt.title('ln(I) gegen $1 / (k_B T)$')
plt.grid(True)
plt.legend()
plt.show()
plt.clf()

#Auswahl der Daten
fit_inv_T1 = np.linspace(2.8e+20,3.1e+20 , 100)
fit_ln_I1 = np.interp(fit_inv_T1, inv_T1, ln_I1)

#Lineare Regression für Reihe 1
fit_params_w1, fit_cov_w1 = np.polyfit(fit_inv_T1, fit_ln_I1, 1, cov=True)
alpha_w1 = ufloat(fit_params_w1[0], np.sqrt(fit_cov_w1[0][0]))
beta_w1 = ufloat(fit_params_w1[1], np.sqrt(fit_cov_w1[1][1]))

#Aktiwierungsenergie für Reihe 1
a1 = -alpha_w1 / k_B
b1 = beta_w1
W1_eV = a1 /  constants.e  # Umrechnung von Joule in eV
print(f'a1 = {a1:.2e} J')
print(f'b1 = {b1:.2e}')
print(f'Aktivierungsenergie W1 = {W1_eV:.2e} eV')

# Plot für ln(I1) gegen 1/(kB*T1) mit Ausgleichsgerade 
plt.scatter(inv_T1, ln_I1, color='purple', marker='x', label='ln(I_1) gegen 1 / (k_B T)')
fit_line_y = np.polyval(fit_params_w1, fit_inv_T1)
plt.plot(fit_inv_T1, fit_line_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / (k_B T) 1/J$')
plt.ylabel(r'$ln(I)$ (Einheitenlos)')
plt.title('ln(I) gegen $1 / (k_B T)$')
plt.grid(True)
plt.legend()
plt.show()
plt.clf()