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
z_nom, Magnetfeld_nom = np.genfromtxt('Magnetfeld.csv', delimiter=',', skip_header=1, unpack='True')
Magnetfeldfehler = np.ones(len(Magnetfeld_nom))*0.5
Magnetfeld= unp.uarray(Magnetfeld_nom, Magnetfeldfehler)

# Plot the results
#plt.figure(figsize=(10, 6))
plt.errorbar(z_nom, Magnetfeld_nom, yerr=Magnetfeldfehler, fmt= "*", color="purple", label= "Messwerte")
plt.xlabel('z/mm')
plt.ylabel('B/mT')
plt.title('Magnetfeldstärke')
plt.legend()
plt.grid(True)
plt.show()

# Save the plot as an image file 
#plt.savefig('/Users/celinawieberg/Documents/Praktikum/FP/v46/plot.png')
plt.clf()
Probe1Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe1Winkel1.csv', delimiter=',', skip_header=1)
Probe1Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe1Winkel2.csv', delimiter=',', skip_header=1)
Probe2Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe2Winkel1.csv', delimiter=',', skip_header=1)
Probe2Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe2Winkel2.csv', delimiter=',', skip_header=1)
Probe3Winkel1 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe3Winkel1.csv', delimiter=',', skip_header=1)
Probe3Winkel2 = np.genfromtxt('/Users/celinawieberg/Documents/Praktikum/FP/v46/Probe3Winkel2.csv', delimiter=',', skip_header=1)

Wellenlaenge = Probe1Winkel1[:, 1]

theta1_1 = Probe1Winkel1[:, 0] * 0.0174533 #Translate deg to rad
thetafehler= np.ones(len(theta1_1))*0.005*0.0174533

theta1_1fehler=unp.uarray(theta1_1,thetafehler) 
theta2_1 = Probe1Winkel2[:, 0] * 0.0174533
theta2_1fehler=unp.uarray(theta2_1,thetafehler) 
thetafrei1= (theta1_1fehler-theta2_1fehler)/(2)
thetafrei1 = unp.uarray(np.abs(unp.nominal_values(thetafrei1)),unp.std_devs(thetafrei1))



theta1_2 = Probe2Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta1_2fehler=unp.uarray(theta1_2,thetafehler) 
theta2_2 = Probe2Winkel2[:, 0] * 0.0174533
theta2_2fehler=unp.uarray(theta2_2,thetafehler) 
thetafrei2= (theta1_2fehler-theta2_2fehler)/(2)
thetafrei2 = unp.uarray(np.abs(unp.nominal_values(thetafrei2)),unp.std_devs(thetafrei2))

theta1_3 = Probe3Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta1_3fehler=unp.uarray(theta1_3,thetafehler) 
theta2_3 = Probe3Winkel2[:, 0] * 0.0174533
theta2_3fehler=unp.uarray(theta2_3,thetafehler) 
thetafrei3= (theta1_3fehler-theta2_3fehler)/(2)
thetafrei3 = unp.uarray(np.abs(unp.nominal_values(thetafrei3)),unp.std_devs(thetafrei3))





plt.errorbar(Wellenlaenge, unp.nominal_values(thetafrei1), yerr=unp.std_devs(thetafrei1), fmt="*", color="purple", label="Probe 1")
plt.errorbar(Wellenlaenge, unp.nominal_values(thetafrei2), yerr=unp.std_devs(thetafrei2), fmt="*", color="green", label="Probe 2")
plt.errorbar(Wellenlaenge, unp.nominal_values(thetafrei3), yerr=unp.std_devs(thetafrei3), fmt="*", color="blue", label="Probe 3")
plt.xlabel(r'$\lambda / \mu$m')
plt.ylabel(r'$\theta$ /rad/m')
plt.legend(loc="best")
plt.grid(True)
plt.show()
plt.clf()

#differenzen der dotierten und undotierten probe

thetadiff1 = (thetafrei1-thetafrei3)
thetadiff2 = (thetafrei2-thetafrei3)

plt.errorbar(Wellenlaenge, unp.nominal_values(thetadiff1),yerr=unp.std_devs(thetadiff1), fmt= "*", color ="purple", label = "Farady Rotation der ersten Probe")
plt.errorbar(Wellenlaenge, unp.nominal_values(thetadiff2),yerr=unp.std_devs(thetadiff2),fmt= "*", color= "green", label = "Faraday Rotation der zweiten Probe")
plt.xlabel(r'$\lambda / \mu$m')
plt.ylabel(r'$\theta_{dot}-\theta_{undot}$/m') 
plt.legend(loc="best")
plt.grid(True)
plt.show()
plt.clf()



from scipy.stats import linregress
# Square of the Wellenlaenge 
Wellenlaenge_squared = Wellenlaenge ** 2 
#Linear Regression for thetadiff1 
slope1, intercept1, r_value1, p_value1, std_err1 = linregress(Wellenlaenge_squared, unp.nominal_values(thetadiff1)) 

# Linear Regression for thetadiff2 
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(Wellenlaenge_squared, unp.nominal_values(thetadiff2)) 

# Extend the regression lines over the entire range of data 
x_fit = np.linspace(min(Wellenlaenge_squared), max(Wellenlaenge_squared), 200)
line1 = slope1 * x_fit + intercept1 
line2 = slope2 * x_fit + intercept2

plt.errorbar(Wellenlaenge_squared ,unp.nominal_values(thetadiff1),yerr=unp.std_devs(thetadiff1), fmt= "*", color ="purple", label = "Farady Rotation der ersten Probe")
plt.errorbar(Wellenlaenge_squared, unp.nominal_values(thetadiff2),yerr=unp.std_devs(thetadiff2),fmt= "*", color= "green", label = "Faraday Rotation der zweiten Probe")
plt.plot(x_fit,line1, "-", color="purple") 
plt.plot(x_fit,line2, "-", color="green")
plt.xlabel(r'$(\lambda)^2 / (\mu m^2)$')
plt.ylabel(r'$\theta_{dot}-\theta_{undot}$/m') 
plt.legend(loc="best")
plt.grid(True)
plt.show()
plt.clf()


# Steigung (Slope) ist m

m_1 = slope
print(f"Proportionalitätsfaktor m_1: {m_1}")
print(std_err)

#Lineare Regression Probe 1:
slope, intercept, r_value, p_value, std_err = linregress(Wellenlaenge_squared ,thetadiff2)

m_1 = slope1
m1 = ufloat(slope1, std_err1)
print(f"Proportionalitätsfaktor m_1: {m1}") 

K1 = 2.39977*10**-63
m_eff1 = unp.sqrt(K1 / m1)
print(f"Effektive Masse m_eff1: {m_eff1}")


# Steigung (Slope) ist m
m_2 = slope2
m2=ufloat(slope2, std_err2)
print(f"Proportionalitätsfaktor m_2: {m2}")
K2=1.9598199*10**-62
m_eff2 = unp.sqrt(K2 /m2)
print(f"Effektive Masse m_eff2: {m_eff2}")