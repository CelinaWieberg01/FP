import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
plt.switch_backend('TkAgg')

#Datensatz 1
t1, T1 , I1= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten1.txt", unpack=True)
I1= I1*10 #umrechnung auf pA
T1=T1+273.15 #umrechnung auf Kelvin

#Linearer Fit für den ersten Untergrund
fit_T1 = np.concatenate((T1[12:17], T1[-16:-9]))  # Anfangs- und Endwerte für den Fit
fit_I1 = np.concatenate((I1[12:17], I1[-16:-9]))
fit_params1, fit_cov1 = np.polyfit(fit_T1, fit_I1, 1, cov=True)
alpha1 = ufloat(fit_params1[0], np.sqrt(fit_cov1[0][0]))
beta1 = ufloat(fit_params1[1], np.sqrt(fit_cov1[1][1]))
fit_line1_x = np.linspace(min(T1), max(T1), 100)
fit_line1_y = np.polyval(fit_params1, fit_line1_x)

#Original Plot für die erste Messreihe
plt.scatter(T1, I1, color='purple', marker='x', label='Erste Heizrate')
plt.plot(fit_line1_x, fit_line1_y, color='green', linestyle='-', label='Untergrund Fit') 
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Der Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 300)
plt.ylim(0, 5)
plt.show()
plt.clf()

#Datensatz 2
t2, T2 , I2= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten2.txt", unpack=True)
I2= I2*10 # umrechnung auf pA
T2=T2+273.15 # umrechnung auf Kelvin

#Linearer Fit für den zweiten Untergrund
fit_T2 = np.concatenate((T2[12:17], T2[-16:-9]))  # Anfangs- und Endwerte für den Fit
fit_I2 = np.concatenate((I2[12:17], I2[-16:-9]))
fit_params2, fit_cov2= np.polyfit(fit_T2, fit_I2, 1, cov=True)
alpha2 = ufloat(fit_params2[0], np.sqrt(fit_cov2[0][0]))
beta2 = ufloat(fit_params2[1], np.sqrt(fit_cov2[1][1]))
fit_line2_x = np.linspace(min(T2), max(T2), 100)
fit_line2_y = np.polyval(fit_params2, fit_line2_x)


# Linearer Fit für die zweite Heizrate
heizrate2_parasm= np.polyfit(t2, T2, 1)  
heizrate2 = heizrate2_parasm[0]

#Original Plot für die zweite Messreihe
plt.scatter(T2, I2, color='purple', marker='x', label='Erste Heizrate') 
plt.plot(fit_line2_x, fit_line2_y, color='green', linestyle='-', label='Untergrund Fit')
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Der Depolarisationsstrom für die zweite Heizrate')
plt.grid(True)
plt.xlim(200, 300)
plt.ylim(0, 5)
plt.show()
plt.clf()

#Angabe Heizraten
print(f'Heizrate für die erste Messung: {alpha1.nominal_value:.2e} K/min')
print(f'Heizrate für die zweite Messung: {alpha2.nominal_value:.2e} K/min')

#Datenbereinigung Reihe 1
untergrund1 = alpha1.nominal_value * T1 + beta1.nominal_value
bereinigt1 = I1-untergrund1
print(f'Erste Messung: alpha1 = {alpha1:.2e} pA/K, beta1 = {beta1:.2e} A')

#Datenbereinigung Reihe 2
untergrund2 = alpha2.nominal_value * T2 + beta2.nominal_value
bereinigt2 = I2-untergrund2
print(f'Zweite Messung: alpha2 = {alpha2:.2e} pA/K, beta2 = {beta2:.2e} A')

#Plot der bereinigten Daten für Reihe 1
plt.scatter(T1, bereinigt1, color='purple', marker='x', label='Bereinigter erster Datensatz durch abziehen des Untergrunds')
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Bereinigter Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 340)
plt.ylim(0, max(bereinigt1))
plt.show()
plt.clf()

#Plot der bereinigten Daten für Reihe 2
plt.scatter(T1, bereinigt1, color='purple', marker='x', label='Bereinigter zweiter Datensatz durch abziehen des Untergrunds')
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Bereinigter Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 340)
plt.ylim(0, max(bereinigt2))
plt.show()
plt.clf()

#Aktivierungsenergie über die Stromdichte
