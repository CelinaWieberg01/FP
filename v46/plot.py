import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('TkAgg')  # oder ein anderes GUI-Backend wie 'Qt5Agg' oder 'Agg'
import matplotlib.pyplot as plt


# Read the CSV file
#data = np.genfromtxt('daten.csv', delimiter=',', skip_header=1)
Magnetfeld = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Magnetfeld.csv', delimiter=',', skip_header=1)


# Extract columns
x_wert = Magnetfeld[:, 0]
y_wert = Magnetfeld[:, 1]


# Plot the results
#plt.figure(figsize=(10, 6))
plt.plot(x_wert, y_wert, "*", color="purple")
plt.xlabel('z/mm')
plt.ylabel('B/mT')
plt.title('Magnetfeldstärke')
plt.legend()
plt.grid(True)
#plt.show()2

# Save the plot as an image file 
#plt.savefig('/Users/celinawieberg/Documents/Praktikum/FP/v46/plot.png')
plt.clf()
Probe1Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe1Winkel1.csv', delimiter=',', skip_header=1)
Probe1Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe1Winkel2.csv', delimiter=',', skip_header=1)
Probe2Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe2Winkel1.csv', delimiter=',', skip_header=1)
Probe2Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe2Winkel2.csv', delimiter=',', skip_header=1)
Probe3Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe3Winkel1.csv', delimiter=',', skip_header=1)
Probe3Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe3Winkel2.csv', delimiter=',', skip_header=1)


theta1_1 = Probe1Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta2_1 = Probe1Winkel2[:, 0] * 0.0174533
Wellenlaenge = Probe1Winkel1[:, 1]
thetafrei1=(theta1_1-theta2_1)/(1.36 *2)

theta1_2 = Probe2Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta2_2 = Probe2Winkel2[:, 0] * 0.0174533
thetafrei2=(theta1_2-theta2_2)/(1.296 *2)

theta1_3 = Probe3Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta2_3 = Probe3Winkel2[:, 0] * 0.0174533
thetafrei3=(theta1_3-theta2_3)/(5.11 *2)





plt.plot(thetafrei1, "*", color="purple")
plt.plot(thetafrei2, "*", color="green")
plt.plot(thetafrei3, "*", color="blue")
plt.xlabel(r'Wellenlänge /$\mu$m')
plt.ylabel(r'$\theta$ /rad')
plt.title('Magnetfeldstärke')
plt.legend()
plt.grid(True)
plt.show()
