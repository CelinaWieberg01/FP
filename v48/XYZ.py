import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

def linfit(x, m, b):
    return m*x+b

t1, T1 , I1= np.genfromtxt("Daten1.txt", unpack=True)
t2, T2 , I2= np.genfromtxt("Daten2.txt", unpack=True)

I1= I1*10
T1=T1+273.15

I2= I2*10
T2=T2+273.15

# Heizraten
plt.scatter(t1, T1, label="Messreihe 1")
plt.scatter(t2, T2, label="Messreihe 2")

params1,  cov1 = curve_fit(linfit,  t1, T1)
params2, cov2 = curve_fit(linfit, t2, T2)

xx = np.linspace(min(t1[0], t2[0]), max(t1[-1], t2[-1]), 1000)
y1 = linfit(xx, *params1)
y2 = linfit(xx, *params2)

plt.plot(xx, y1, label = "Heizrate 1")
plt.plot(xx, y2, label="Heizrate 2")
print("Heizrate 1 über Fit = ", params1[0], " +- ", np.diag(np.sqrt(cov1))[0])
print("Heizrate 2 über Fit = ", params2[0], " +- ", np.diag(np.sqrt(cov2))[0])
plt.legend()
plt.grid("on")
plt.figure()

"""
Heizrate 1 über Fit =  1.2420890575344499  +-  0.019833746560421136
Heizrate 2 über Fit =  1.6799985616114808  +-  0.019616136972691517
"""


#Linearer Fit für den ersten Untergrund
fit_T1 = np.concatenate((T1[:10], T1[-16:]))  # Anfangs- und Endwerte für den Fit
fit_I1 = np.concatenate((I1[:10], I1[-16:]))

params_u1, cov_u1 = curve_fit(linfit, fit_T1, fit_I1)
xx_fit1 = np.linspace(fit_T1[0], fit_T1[-1], 1000)
yy_fit1 = linfit(xx_fit1, *params_u1)

print("Parameter für Untergrund = ", params_u1, " +- ", np.sqrt(np.diag(cov_u1)))


# Linearer Fit für die erste Heizrate
heizzrate1_params1 = np.polyfit(t1, T1, 1)  # Grad 1 für linearen Fit
heizrate1 = heizzrate1_params1[0]


plt.scatter(T1, I1, color='purple', marker='x', label='Erste Heizrate')
plt.plot(xx_fit1, yy_fit1, color='green', linestyle='-', label='Untergrund Fit') 
plt.xlabel('T in K')
plt.ylabel('I in pA')
plt.title('Der Depolarisationsstrom für die erste Heizrate')
plt.grid(True)
plt.legend()
plt.figure()


#Linearer Fit für den zweiten Untergrund
fit_T2 = np.concatenate((T2[:12], T2[-9:]))  # Anfangs- und Endwerte für den Fit
fit_I2 = np.concatenate((I2[:12], I2[-9:]))

params_u2, cov_u2 = curve_fit(linfit, fit_T2, fit_I2)
xx_fit2 = np.linspace(fit_T2[0], fit_T2[-1], 1000)
yy_fit2 = linfit(xx_fit2, *params_u2)

print("Parameter für Untergrund = ", params_u2, " +- ", np.sqrt(np.diag(cov_u2)))


# Linearer Fit für die zweite Heizrate
heizrate2_parasm= np.polyfit(t2, T2, 1)  
heizrate2 = heizrate2_parasm[0]


plt.scatter(T2, I2, color='purple', marker='x', label='Erste Heizrate') 
plt.plot(xx_fit2, yy_fit2, color='green', linestyle='-', label='Untergrund Fit')
plt.xlabel('T in K')
plt.ylabel('I in pA')
plt.title('Der Depolarisationsstrom für die zweite Heizrate')
plt.grid(True)
plt.legend()
plt.figure()

# Anzeigen der Heizraten
print(f'Heizrate für die erste Messung: {heizrate1:.2f} K/min')
print(f'Heizrate für die zweite Messung: {heizrate2:.2f} K/min')


"""
Von Celina
Fit Untergrund 1 =  [ 0.0122142  -2.31501693]
Fit Untergrund 2 =  [ 0.02114334 -4.35632341]

von mir
Parameter für Untergrund =  [ 0.0122142  -2.31501693]  +-  [0.00037518 0.09826862]
Parameter für Untergrund =  [ 0.02114334 -4.35632341]  +-  [0.00048343 0.12130725]


"""