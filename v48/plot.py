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
#plt.show()
plt.clf()

#Angabe Heizrate 2
print(f'Heizrate für die zweite Messung: {alpha2.nominal_value:.2e} K/min')

#Datenbereinigung Reihe 2
untergrund2 = alpha2.nominal_value * T2 + beta2.nominal_value
bereinigt2 = I2-untergrund2
print(f'Zweite Messung: alpha2 = {alpha2:.2e} pA/K, beta2 = {beta2:.2e} A')

#Plot der bereinigten Daten für Reihe 2
plt.scatter(T1, bereinigt1, color='purple', marker='x', label='Bereinigter zweiter Datensatz durch abziehen des Untergrunds')
plt.xlabel(r'$T$ in K')
plt.ylabel(r'$I$ in pA')
plt.title('Bereinigter Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.xlim(200, 340)
plt.ylim(0, max(bereinigt2))
#plt.show()
plt.clf()

#Aktivierungsenergie über den Polarisationsansatz
k_B = constants.Boltzmann

#erechnung des Logarithmus von I und der Inversen Temperatur für Messreihe 1
I1= I1 * 1e-11#auf Ampere
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



"""
# Plot für ln(I1) gegen 1/(kB*T1) für den kleinen Temperaturbereich
plt.scatter(inv_T1, ln_I1, color='purple', marker='x', label='ln(I_1) gegen 1 / (k_B T)')
fit_line_y = np.polyval(fit_params_W1, fit_line_x)
plt.plot(fit_line_x, fit_line_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / (k_B T) \cdot 10^{20} \, \text{J}^{-1}$')
plt.ylabel(r'$ln(I)$ (Einheitenlos)')
plt.title('ln(I) gegen $1 / (k_B T) \cdot 10^{20} \, \text{J}^{-1}$ für die erste Heizrate')
plt.grid(True)
plt.legend()
plt.show()
plt.clf()"""

"""# Plot für ln(I1) gegen 1/ K_b T
plt.scatter(1 / (k_B * T1) , np.log(I1), color='red', marker='x', label='ln(I) gegen 1 / kBT')
fit_line_W1_x = np.linspace(min(inv_T1), max(inv_T1), 100)
fit_line_W1_y = np.polyval(fit_params_W1, fit_line_W1_x)
plt.plot(fit_line_W1_x, fit_line_W1_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / k_b T$ in 1/J')
plt.ylabel(r'$ln(I)$ in 1/A')
plt.title('ln(I) gegen $1 / k_b T$ 1/Jfür die erste Heizrate')
plt.grid(True)
#plt.legend()
plt.show()
plt.clf()"""

"""# Messreihe 2
ln_I2 = np.log(I2[14:-19])  # Logarithmus des Stroms
inv_T2 = 1 / (k_B * T2[14:-19])  # Inverse Temperatur
print(ln_I2)
fit_params_W2, fit_cov_W2 = np.polyfit(inv_T2, ln_I2, 1, cov=True)  # Lineare Regression
alpha_W2 = ufloat(fit_params_W2[0], np.sqrt(fit_cov_W2[0][0]))
beta_W2 = ufloat(fit_params_W2[1], np.sqrt(fit_cov_W2[1][1]))

W2 = -alpha_W2 * k_B 

print(f'Aktivierungsenergie W2 = {W2:.2f} eV')

# Plot für ln(I2) gegen 10^-18 J / K_b T
plt.scatter(1 / (k_B * T2) , np.log(I2), color='red', marker='x', label='ln(I_2) gegen 1 / k_b T')
fit_line_W2_x = np.linspace(44.4 ,51, 100)
fit_line_W2_y = np.polyval(fit_params_W2, fit_line_W2_x)
plt.plot(fit_line_W2_x, fit_line_W2_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / k_b T$ eV')
plt.ylabel(r'$ln(I)$ in pA')
plt.title('ln(I) gegen $1 / k_b$')
plt.grid(True)
plt.show()
plt.clf()"""

"""# Aktivierungsenergie und Relaxationszeit über den Ansatz der Stromdichte
def grafik_integration(I, T):
    return np.cumsum(I * np.gradient(T))

# Hilfsfunktion zur Berechnung von f(T)
def berechne_f_T(I, T):
    P = grafik_integration(I, T)
    f_T = np.log(P) - np.log(I)
    return f_T

# Messreihe 1
f_T1 = berechne_f_T(I1, T1)
print("Werte von ln(f(T)) für Messreihe 1:", f_T1)

# Plot ohne Regression anzeigen
plt.scatter(1 / (k_B * T1) , f_T1, color='purple', marker='x', label='ln(f(T_1)) gegen 1 / k_b T')
plt.xlabel(r'$1 / k_b T$ [eV]')
plt.ylabel(r'$ln(f(T))$')
plt.title('ln(f(T)) gegen $1 / k_b T$ [eV] für die erste Heizrate')
plt.grid(True)
plt.legend()
plt.show()
plt.clf()

# Definition des Bereichs für die Regression für Messreihe 1
fit_line_tau1_x = np.linspace(40, 60, 100)  # Manuelle Definition des Bereichs
relevante_indices1 = np.where((1 / (k_B * T1) >= 0) & (1 / (k_B * T1) <= 100))
print("ggggggggggggggg",relevante_indices1)
inv_T1_relevant = 1 / (k_B *T1[relevante_indices1]) 
f_T1_relevant = f_T1[relevante_indices1]

print("Relevante Werte für Messreihe 1 - inv_T1_relevant:", inv_T1_relevant)
print("Relevante Werte für Messreihe 1 - f_T1_relevant:", f_T1_relevant)

fit_params_tau1, fit_cov_tau1 = np.polyfit(inv_T1_relevant, f_T1_relevant, 1, cov=True)  # Lineare Regression
fit_line_tau1_y = np.polyval(fit_params_tau1, fit_line_tau1_x)

# Plot der Regression
plt.scatter(1 / (k_B * T1) , f_T1, color='red', marker='x', label='ln(f(T_1)) gegen 1 / k_b T')
plt.plot(fit_line_tau1_x, fit_line_tau1_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / k_b T$ [eV]')
plt.ylabel(r'$ln(f(T))$')
plt.title('ln(f(T)) gegen $1 / k_b T$ [eV] für die erste Heizrate')
plt.grid(True)
plt.legend()
plt.xlim(40, 60)
#plt.ylim(-2, 6)
plt.show()
plt.clf()

# Messreihe 2
f_T2 = berechne_f_T(I2, T2)
print("Werte von ln(f(T)) für Messreihe 2:", f_T2)

# Plot ohne Regression anzeigen
plt.scatter(1 / (k_B * T2) , f_T2, color='purple', marker='x', label='ln(f(T_2)) gegen 1 / k_b T')
plt.xlabel(r'$1 / k_b T$ [eV]')
plt.ylabel(r'$ln(f(T))$')
plt.title('ln(f(T)) gegen $1 / k_b T$ [eV] für die zweite Heizrate')
plt.grid(True)
plt.legend()
#plt.show()
plt.clf()

# Definition des Bereichs für die Regression für Messreihe 2
fit_line_tau2_x = np.linspace(46.5, 49.5, 100)  # Manuelle Definition des Bereichs
relevante_indices2 = np.where((1 / (k_B * T2)>= 46.5) & (1 / (k_B * T2) <= 49.5))
inv_T2_relevant = 1 / T2[relevante_indices2]
f_T2_relevant = f_T2[relevante_indices2]

print("Relevante Werte für Messreihe 2 - inv_T2_relevant:", inv_T2_relevant)
print("Relevante Werte für Messreihe 2 - f_T2_relevant:", f_T2_relevant)

fit_params_tau2, fit_cov_tau2 = np.polyfit(inv_T2_relevant, f_T2_relevant, 1, cov=True)  # Lineare Regression
fit_line_tau2_y = np.polyval(fit_params_tau2, fit_line_tau2_x)

# Plot der Regression
plt.scatter(1 / (k_B * T2) , f_T2, color='red', marker='x', label='ln(f(T_2)) gegen 1 / k_b T')
plt.plot(fit_line_tau2_x, fit_line_tau2_y, color='green', linestyle='-', label='Fit')
plt.xlabel(r'$1 / k_b T$ [eV]')
plt.ylabel(r'$ln(f(T))$')
plt.title('ln(f(T)) gegen $1 / k_b T$ [eV] für die zweite Heizrate')
plt.grid(True)
plt.legend()
plt.xlim(46.5, 49.5)
plt.ylim(-2, 6)
#plt.show()
plt.clf()"""