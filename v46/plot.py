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
#plt.show()

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
thetafrei1= np.abs(theta1_1-theta2_1)/(1.36 *2)

theta1_2 = Probe2Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta2_2 = Probe2Winkel2[:, 0] * 0.0174533
thetafrei2= np.abs(theta1_2-theta2_2)/(1.296 *2)

theta1_3 = Probe3Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta2_3 = Probe3Winkel2[:, 0] * 0.0174533
thetafrei3= np.abs(theta1_3-theta2_3)/(5.11 *2)





plt.plot(Wellenlaenge, thetafrei1, "*", color="purple", label="Probe 1")
plt.plot(Wellenlaenge, thetafrei2, "*", color="green", label="Probe 2")
plt.plot(Wellenlaenge, thetafrei3, "*", color="blue", label="Probe 3")
plt.xlabel(r'$\lambda / \mu$m')
plt.ylabel(r'$\theta$ /rad')
plt.legend(loc="best")
plt.grid(True)
plt.show()
plt.clf()

#differenzen der dotierten und undotierten probe

thetadiff1 = np.abs(thetafrei1-thetafrei3)
thetadiff2 = np.abs(thetafrei2-thetafrei3)

plt.plot(Wellenlaenge, thetadiff1, "*", color ="purple", label = "Farady Rotation der ersten Probe")
plt.plot(Wellenlaenge, thetadiff2, "*", color= "green", label = "Faraday Rotation der zweiten Probe")
plt.xlabel(r'$\lambda / \mu$m')
plt.ylabel(r'$\theta_{dot}-\theta_{undot}$') 
plt.legend(loc="best")
plt.grid(True)
plt.show()
plt.clf()



from scipy.stats import linregress
# Lineare Regression für thetadiff1
slope1, intercept1, r_value1, p_value1, std_err1 = linregress(Wellenlaenge*Wellenlaenge, thetadiff1)
line1 = slope1 * np.array(Wellenlaenge*Wellenlaenge) + intercept1

# Lineare Regression für thetadiff2
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(Wellenlaenge*Wellenlaenge, thetadiff2)
line2 = slope2 * np.array(Wellenlaenge*Wellenlaenge) + intercept2

plt.plot(Wellenlaenge*Wellenlaenge, thetadiff1, "*", color ="purple", label = "Farady Rotation der ersten Probe")
plt.plot(Wellenlaenge*Wellenlaenge, thetadiff2, "*", color= "green", label = "Faraday Rotation der zweiten Probe")
plt.plot(line1, "-", color="purple", label=f"Lineare Regression Probe 1 (R²={r_value1**2:.2f})") 
plt.plot(line2, "-", color="green", label=f"Lineare Regression Probe 2 (R²={r_value2**2:.2f})")
plt.xlabel(r'$(\lambda)^2 / (\mu m^2)$')
plt.ylabel(r'$\theta_{dot}-\theta_{undot}$') 
plt.legend(loc="best")
plt.grid(True)
plt.show()
