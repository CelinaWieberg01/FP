import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl



# Read the CSV file
#data = np.genfromtxt('daten.csv', delimiter=',', skip_header=1)
z_nom, Magnetfeld_nom = np.genfromtxt('Magnetfeld.csv', delimiter=',', skip_header=1, unpack='True')
Magnetfeldfehler = np.ones(len(Magnetfeld_nom))*0.5
Magnetfeld= unp.uarray(Magnetfeld_nom, Magnetfeldfehler)

# Plot the results
#plt.figure(figsize=(10, 6))
plt.errorbar(z_nom, Magnetfeld_nom, yerr=Magnetfeldfehler, fmt= "*", color="purple", label= "Messwerte")
plt.xlabel(r'$z$ in $\si{\milli\meter}$')
plt.ylabel(r'$B$ in $\si{\milli\tesla}$')
plt.title('Magnetfeldstärke')
plt.legend()
plt.grid(True)
plt.savefig("plots/magnetfeld.pdf")
plt.figure()

# Save the plot as an image file 
#plt.savefig('/Users/celinawieberg/Documents/Praktikum/FP/v46/plot.png')
plt.clf()
Probe1Winkel1 = np.genfromtxt('Probe1Winkel1.csv', delimiter=',', skip_header=1)
Probe1Winkel2 = np.genfromtxt('Probe1Winkel2.csv', delimiter=',', skip_header=1)
Probe2Winkel1 = np.genfromtxt('Probe2Winkel1.csv', delimiter=',', skip_header=1)
Probe2Winkel2 = np.genfromtxt('Probe2Winkel2.csv', delimiter=',', skip_header=1)
Probe3Winkel1 = np.genfromtxt('Probe3Winkel1.csv', delimiter=',', skip_header=1)
Probe3Winkel2 = np.genfromtxt('Probe3Winkel2.csv', delimiter=',', skip_header=1)

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
plt.xlabel(r'$\lambda$ in $\si{\micro\meter}$')
plt.ylabel(r'$\theta$ in $\si{\radian}$')
plt.legend(loc="best")
plt.grid(True)
plt.savefig("plots/raw_data.pdf")
plt.figure()

#differenzen der dotierten und undotierten probe

thetadiff1 = (thetafrei1-thetafrei3)
thetadiff2 = (thetafrei2-thetafrei3)

plt.errorbar(Wellenlaenge, unp.nominal_values(thetadiff1),yerr=unp.std_devs(thetadiff1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge, unp.nominal_values(thetadiff2),yerr=unp.std_devs(thetadiff2),fmt= "*", color= "green", label = "Probe 2")
plt.xlabel(r'$\lambda$ in $\si{\micro\meter}$')
plt.ylabel(r'$\theta_\text{dot}-\theta_\text{undot}$ in $\si{\radian}$') 
plt.legend(loc="best")
plt.grid(True)
#plt.show()
plt.clf()
plt.figure()



from scipy.stats import linregress
# Square of the Wellenlaenge 
Wellenlaenge_squared = Wellenlaenge ** 2
x_fit = np.linspace(min(Wellenlaenge_squared), max(Wellenlaenge_squared), 200)

#Linear Regression for thetadiff1 
res1 = linregress(Wellenlaenge_squared, unp.nominal_values(thetadiff1)) 
print(f"slope1 = {res1.slope} +- {res1.stderr}, intercept1 = {res1.intercept} +- {res1.intercept_stderr}")
line1 = res1.slope * x_fit + res1.intercept 


# Linear Regression for thetadiff2 
res2 = linregress(Wellenlaenge_squared, unp.nominal_values(thetadiff2)) 
print(f"slope1 = {res2.slope} +- {res2.stderr}, intercept1 = {res2.intercept} +- {res2.intercept_stderr}")
line2 = res2.slope * x_fit + res2.intercept




plt.errorbar(Wellenlaenge_squared ,unp.nominal_values(thetadiff1),yerr=unp.std_devs(thetadiff1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge_squared, unp.nominal_values(thetadiff2),yerr=unp.std_devs(thetadiff2),fmt= "*", color= "green", label = "Probe 2")
plt.plot(x_fit,line1, "-", color="purple", label="Fit Probe 1") 
plt.plot(x_fit,line2, "-", color="green", label="Fit Probe 2")
plt.xlabel(r'$\lambda^2$ in $\si{\square\micro\meter}$')
plt.ylabel(r'$\theta_\text{dot}-\theta_\text{undot}$ in $\si{\radian}$') 

plt.legend(loc="best")
plt.grid(True)

plt.savefig("plots/fits_ohne_n.pdf")
plt.figure()

e0 = constants.elementary_charge
eps = constants.epsilon_0 * 12.4
c = constants.speed_of_light
B = 421 * 10**-3
n = 3.5

# Propfaktor
def faktor_festes_n(N):
    k = np.sqrt(e0**3 * N * B / ( 8*np.pi**2 * eps * c**3 * n))
    return k

#Lineare Regression Probe 1:
m1 = ufloat(res1.slope, res1.stderr)
print(f"Proportionalitätsfaktor m_1: {m1}") 

N1 = 1.2e12
k1 = faktor_festes_n(N1)
print("k1 = ", k1)
m_eff_aaron = k1 / unp.sqrt(m1)
print(f"m_eff_aaron = {m_eff_aaron}")




# Steigung (Slope) ist m
m2=ufloat(res2.slope, res2.stderr)
print(f"Proportionalitätsfaktor m_2: {m2}")

N2 = 2.8e12
k2 = faktor_festes_n(N2)
print("k2 = ", k2)
m_eff2 = k2 / unp.sqrt(m2)
print(f"Effektive Masse m_eff2: {m_eff2}")










# now lets do this again for the different refractive indices, for which we fit a function to wavelengths squared divided by respective refractive indices

# to find refractive indices of wavelengths, we use the sellmeyer equation for gallium arsenide

def sellmeyer(welle_squared):
    n = np.sqrt(8.950 + 2.054/(1 - 0.390 / welle_squared))
    return n

from scipy.stats import linregress
# Square of the Wellenlaenge 
Wellenlaenge_squared = Wellenlaenge**2 
n_array = sellmeyer(Wellenlaenge_squared)

Wellenlaenge_squared_n = Wellenlaenge_squared / n_array

xx = np.linspace(Wellenlaenge_squared_n[0], Wellenlaenge_squared_n[-1], 100)
plt.errorbar(Wellenlaenge_squared_n,unp.nominal_values(thetadiff1),yerr=unp.std_devs(thetadiff1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge_squared_n, unp.nominal_values(thetadiff2),yerr=unp.std_devs(thetadiff2),fmt= "*", color= "green", label = "Probe 2")

#Linear Regression for thetadiff1 
res1n = linregress(Wellenlaenge_squared_n, unp.nominal_values(thetadiff1)) 
m1n = ufloat(res1n.slope, res1n.stderr)
b1n = ufloat(res1n.intercept, res1n.intercept_stderr)
print("m1 mit n beachtet = ", m1n)
print("b1 mit n beachtet = ", b1n)

fit1 = res1n.slope * xx + res1n.intercept
plt.plot(xx, fit1, color="purple", label="Fit Probe 1")

#Linear Regression for thetadiff1 
res2n = linregress(Wellenlaenge_squared_n, unp.nominal_values(thetadiff2)) 
m2n = ufloat(res2n.slope, res2n.stderr)
b2n = ufloat(res2n.intercept, res2n.intercept_stderr)
print("m2 mit n beachtet = ", m2n)
print("b2 mit n beachtet = ", b2n)

fit2 = res2n.slope * xx + res2n.intercept
plt.plot(xx, fit2, color="green", label="Fit Probe 2")


plt.xlabel(r'$\frac{\lambda^2}{n}$ in $\si{\square\micro\meter}$')
plt.ylabel(r'$\theta_\text{dot}-\theta_\text{undot}$ in $\si{\radian}$') 

plt.legend(loc="best")
plt.grid(True)

plt.savefig("plots/fits_mit_n.pdf")
plt.figure()

#Lineare Regression Probe 1:
N1 = 1.2e12
k1 = faktor_festes_n(N1)
print("k1 = ", k1)
m_eff_aaron = k1 / unp.sqrt(m1n)
print(f"m_eff_aaron = {m_eff_aaron}")




# Steigung (Slope) ist m
N2 = 2.8e12
k2 = faktor_festes_n(N2)
print("k2 = ", k2)
m_eff2 = k2 / unp.sqrt(m2n)
print(f"Effektive Masse m_eff2: {m_eff2}")



"""
slope1 = 0.02958028375759921 +- 0.01557348533312387, intercept1 = -0.1599743298610714 +- 0.0679915882816386
k1 =  5.04141107676392e-32
m_eff_aaron = (2.9+/-0.8)e-31

slope1 = 0.03865883949704651 +- 0.015621923091931628, intercept1 = -0.16237578851477455 +- 0.06820306054258059
k2 =  7.700882622886497e-32
Effektive Masse m_eff2: (3.9+/-0.8)e-31


n_array =  [3.4779243  3.41067826 3.38700048 3.36391003 3.3520268  3.34546196
 3.34132547 3.33759507 3.33538273]

m1 mit n beachtet =  0.10+/-0.05
b1 mit n beachtet =  -0.16+/-0.07
k1 =  5.04141107676392e-32
m_eff_aaron = (1.6+/-0.4)e-31


m2 mit n beachtet =  0.13+/-0.05
b2 mit n beachtet =  -0.16+/-0.07
k2 =  7.700882622886497e-32
Effektive Masse m_eff2: (2.2+/-0.4)e-31
"""

me = constants.electron_mass
m_lit = 0.067*me
m_ex = np.array([2.9e-31, 3.9e-31, 1.6e-31, 2.2e-31])

for m in m_ex:
    print(abs(m - me)/me)