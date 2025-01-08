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
theta1= (theta1_1fehler-theta2_1fehler)/(2)
plt.errorbar(Wellenlaenge, unp.nominal_values(theta1), yerr=unp.std_devs(theta1), fmt="*", color="purple", label="Probe 1")
L1 = 1.36*10**(-3)
theta1 = unp.uarray(np.abs(unp.nominal_values(theta1)),unp.std_devs(theta1))/L1

theta1_2 = Probe2Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta1_2fehler=unp.uarray(theta1_2,thetafehler) 
theta2_2 = Probe2Winkel2[:, 0] * 0.0174533
theta2_2fehler=unp.uarray(theta2_2,thetafehler) 
theta2= (theta1_2fehler-theta2_2fehler)/(2)
plt.errorbar(Wellenlaenge, unp.nominal_values(theta2), yerr=unp.std_devs(theta2), fmt="*", color="green", label="Probe 2")
L2 = 1.296*10**(-3)
theta2 = unp.uarray(np.abs(unp.nominal_values(theta2)),unp.std_devs(theta2))/L2

theta1_3 = Probe3Winkel1[:, 0] * 0.0174533 #Translate deg to rad
theta1_3fehler=unp.uarray(theta1_3,thetafehler) 
theta2_3 = Probe3Winkel2[:, 0] * 0.0174533
theta2_3fehler=unp.uarray(theta2_3,thetafehler) 
theta3= (theta1_3fehler-theta2_3fehler)/(2)
plt.errorbar(Wellenlaenge, unp.nominal_values(theta3), yerr=unp.std_devs(theta3), fmt="*", color="blue", label="Probe 3")
L3 = 5.11*10**(-3)
theta3 = unp.uarray(np.abs(unp.nominal_values(theta3)),unp.std_devs(theta3))/L3


plt.xlabel(r'$\lambda$ in $\si{\micro\meter}$')
plt.ylabel(r'$\theta$ in $\si{\radian}$')
plt.legend(loc="best")
plt.grid(True)
plt.savefig("plots/raw_data.pdf")
plt.figure()

#differenzen der dotierten und undotierten probe

thetafrei1 = theta1 - theta3
thetafrei2 = theta2 - theta3

plt.errorbar(Wellenlaenge, unp.nominal_values(thetafrei1),yerr=unp.std_devs(thetafrei1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge, unp.nominal_values(thetafrei2),yerr=unp.std_devs(thetafrei2),fmt= "*", color= "green", label = "Probe 2")
plt.xlabel(r'$\lambda$ in $\si{\micro\meter}$')
plt.ylabel(r'$\theta_\text{frei} in \si{radian}$') 
plt.legend(loc="best")
plt.grid(True)
#plt.show()
plt.clf()
plt.figure()


from scipy.stats import linregress

def fit(x, steigung):
    return steigung*x

# Square of the Wellenlaenge 
Wellenlaenge_squared = Wellenlaenge ** 2
x_fit = np.linspace(0, max(Wellenlaenge_squared)+0.25, 200)

#Linear Regression for thetadiff1 
params1, cov1 = curve_fit(fit, Wellenlaenge_squared, unp.nominal_values(thetafrei1), absolute_sigma=True, sigma=unp.std_devs(thetafrei1)) 
print(f"m1 = {params1} +- {np.sqrt(cov1)}")
line1 = params1 * x_fit


# Linear Regression for thetadiff2 
params2, cov2 = curve_fit(fit, Wellenlaenge_squared, unp.nominal_values(thetafrei2), absolute_sigma=True, sigma=unp.std_devs(thetafrei2)) 
print(f"m2 = {params2} +- {np.sqrt(cov2)}")
line2 = params2 * x_fit




plt.errorbar(Wellenlaenge_squared ,unp.nominal_values(thetafrei1),yerr=unp.std_devs(thetafrei1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge_squared, unp.nominal_values(thetafrei2),yerr=unp.std_devs(thetafrei2),fmt= "*", color= "green", label = "Probe 2")
plt.plot(x_fit,line1, "-", color="purple", label="Fit Probe 1") 
plt.plot(x_fit,line2, "-", color="green", label="Fit Probe 2")
plt.xlabel(r'$\lambda^2$ in $10^{-12} \cdot \si{\square\meter}$')
plt.ylabel(r'$\theta_\text{frei, dot}-\theta_\text{frei, undot}$ in $\si{\radian\per\meter}$') 

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
m1 = ufloat(params1, np.sqrt(cov1))

N1 = 1.2e12
k1 = faktor_festes_n(N1)
m_eff1 = k1 / unp.sqrt(m1)
print(f"m_eff1 = {m_eff1}")




# Steigung (Slope) ist m
m2=ufloat(params2, np.sqrt(cov2))

N2 = 2.8e12
k2 = faktor_festes_n(N2)
m_eff2 = k2 / unp.sqrt(m2)
print(f"m_eff2: {m_eff2}")










# now lets do this again for the different refractive indices, for which we fit a function to wavelengths squared divided by respective refractive indices

# to find refractive indices of wavelengths, we use the sellmeyer equation for gallium arsenide

def sellmeyer(welle_squared):
    n = np.sqrt(8.950 + 2.054/(1 - 0.390 / welle_squared))
    return n

def fit(x, steigung):
    return steigung*x

from scipy.stats import linregress
# Square of the Wellenlaenge 
Wellenlaenge_squared = Wellenlaenge**2 
n_array = sellmeyer(Wellenlaenge_squared)

Wellenlaenge_squared_n = Wellenlaenge_squared / n_array

xx = np.linspace(0, Wellenlaenge_squared_n[-1]+0.25, 100)
plt.errorbar(Wellenlaenge_squared_n,unp.nominal_values(thetafrei1),yerr=unp.std_devs(thetafrei1), fmt= "*", color ="purple", label = "Probe 1")
plt.errorbar(Wellenlaenge_squared_n, unp.nominal_values(thetafrei2),yerr=unp.std_devs(thetafrei2),fmt= "*", color= "green", label = "Probe 2")

#Linear Regression for thetadiff1 
params_n1, cov_n1 = curve_fit(fit, Wellenlaenge_squared_n, unp.nominal_values(thetafrei1), absolute_sigma=True, sigma=unp.std_devs(thetafrei1)) 
m_n1 = ufloat(params_n1, np.sqrt(cov_n1))
print("m_n1 = ", m_n1)


fit1 = params_n1*xx
plt.plot(xx, fit1, color="purple", label="Fit Probe 1")

#Linear Regression for thetadiff1 
params_n2, cov_n2 = curve_fit(fit, Wellenlaenge_squared_n, unp.nominal_values(thetafrei2), absolute_sigma=True, sigma=unp.std_devs(thetafrei2)) 
m_n2 = ufloat(params_n2, np.sqrt(cov_n2))
print("m_n2 = ", m_n2)

fit2 = params_n2*xx
plt.plot(xx, fit2, color="green", label="Fit Probe 2")


plt.xlabel(r'$\frac{\lambda^2}{n}$ in $10^{-12} \cdot \si{\square\meter}$')
plt.ylabel(r'$\theta_\text{frei, dot}-\theta_\text{frei, undot}$ in $\si{\radian\per\meter}$') 

plt.legend(loc="best")
plt.grid(True)

plt.savefig("plots/fits_mit_n.pdf")
plt.figure()

#Lineare Regression Probe 1:
N1 = 1.2e12
k1 = faktor_festes_n(N1)
m_n_eff1 = k1 / unp.sqrt(m_n1)
print(f"m_n_eff1 = {m_n_eff1}")




# Steigung (Slope) ist m
N2 = 2.8e12
k2 = faktor_festes_n(N2)
m_n_eff2 = k2 / unp.sqrt(m_n2)
print(f"m_n_eff2: {m_n_eff2}")


me = constants.electron_mass

m_lit = 0.067*me
print("effektive Elektronenmasse = ", m_lit)
m_ex = np.array([m_eff1, m_eff2, m_n_eff1, m_n_eff2])

for m in m_ex:
    print(abs(m - m_lit)/m_lit)


"""
slope1 = 5.480407351580117 +- 3.4385492215926026, intercept1 = 4.755816647218399 +- 15.012209403338295
slope1 = 12.46471975542822 +- 3.5223613173511623, intercept1 = 5.329526407717992 +- 15.3781209116451
Proportionalitätsfaktor m_1: 5.5+/-3.4
k1 =  5.04141107676392e-32
m_eff_aaron = (2.2+/-0.7)e-32
Proportionalitätsfaktor m_2: 12.5+/-3.5
k2 =  7.700882622886497e-32
Effektive Masse m_eff2: (2.18+/-0.31)e-32
m1 mit n beachtet =  18+/-11
b1 mit n beachtet =  5+/-15
m2 mit n beachtet =  41+/-12
b2 mit n beachtet =  6+/-15
k1 =  5.04141107676392e-32
m_eff_aaron = (1.2+/-0.4)e-32
k2 =  7.700882622886497e-32
Effektive Masse m_eff2: (1.20+/-0.17)e-32

abweichungen:
0.6816469593302487
0.5718700487544723
0.8243569430787578
0.758490796733292
"""

### LINEARER FIT OHNE ACHSENABSCHNITT

"""

m1 = [6.45192752] +- [[0.00358479]]
m2 = [13.55343748] +- [[0.00375037]]
m_eff1 = (1.9848+/-0.0006)e-32
m_eff2: (2.09178+/-0.00029)e-32

m_n1 =  21.567+/-0.012
m_n_eff1 = (1.08558+/-0.00030)e-32

m_n2 =  45.316+/-0.013
m_n_eff2: (1.14396+/-0.00016)e-32

effektive Elektronenmasse =  6.103287080005001e-32
0.67481+/-0.00009
0.65727+/-0.00005
0.82213+/-0.00005
0.812566+/-0.000026

"""