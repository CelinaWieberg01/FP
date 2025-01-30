import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
plt.switch_backend('TkAgg')

t1, T1 , I1= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten1.txt", unpack=True)
I1= I1*10
T1=T1+273.15

#Linearer Fit für den ersten Untergrund
fit_T1 = np.concatenate((T1[:10], T1[-16:]))  # Anfangs- und Endwerte für den Fit
fit_I1 = np.concatenate((I1[:10], I1[-16:]))
fit_params1 = np.polyfit(fit_T1, fit_I1, 1)
fit_line1_x = np.linspace(min(T1), max(T1), 100)
fit_line1_y = np.polyval(fit_params1, fit_line1_x)



print(fit_line1_x)
fit_line1_y = np.polyval(fit_params1, fit_line1_x) 

# Linearer Fit für die erste Heizrate
heizzrate1_params1 = np.polyfit(t1, T1, 1)  # Grad 1 für linearen Fit
heizrate1 = heizzrate1_params1[0]


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

t2, T2 , I2= np.genfromtxt("/Users/celinawieberg/Documents/Praktikum/FP/v48/Daten2.txt", unpack=True)
I2= I2*10
T2=T2+273.15
#Linearer Fit für den zweiten Untergrund
fit_T2 = np.concatenate((T2[:12], T2[-9:]))  # Anfangs- und Endwerte für den Fit
fit_I2 = np.concatenate((I2[:12], I2[-9:]))
fit_params2 = np.polyfit(fit_T2, fit_I2, 1)
fit_line2_x = np.linspace(min(T2), max(T2), 100)
fit_line2_y = np.polyval(fit_params2, fit_line2_x)


# Linearer Fit für die zweite Heizrate
heizrate2_parasm= np.polyfit(t2, T2, 1)  
heizrate2 = heizrate2_parasm[0]


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

# Anzeigen der Heizraten
print(f'Heizrate für die erste Messung: {heizrate1:.2f} K/min')
print(f'Heizrate für die zweite Messung: {heizrate2:.2f} K/min')

