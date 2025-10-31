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

I1= I1 * 10**(-11)
T1=T1+273.15

I2= I2 * 10**(-11)
T2=T2+273.15

T1_err = np.ones(len(T1))*0.3
T2_err = np.ones(len(T2))*0.3
I1_err = np.ones(len(I1))*0.05*10**(-11)
I2_err = np.ones(len(I2))*0.05*10**(-11)

T1_true = unp.uarray(T1, T1_err)
T2_true = unp.uarray(T2, T2_err)
I1_true = unp.uarray(I1, I1_err)
I2_true = unp.uarray(I2, I2_err)


# Heizraten
plt.errorbar(t1, T1, yerr=T1_err, color="orange", marker="x", markersize=9, linestyle="none", label="1. Messreihe")
plt.errorbar(t2, T2, yerr=T2_err, color="blue",   marker="+", markersize=9, linestyle="none", label="2. Messreihe")
plt.plot(t1, T1, color="orange", linewidth=3,  alpha=0.4)
plt.plot(t2, T2, color="blue",   linewidth=3,  alpha=0.4)

params1,  cov1 = curve_fit(linfit,  t1, T1)
params2, cov2 = curve_fit(linfit, t2, T2)

xx = np.linspace(min(t1[0], t2[0]), max(t1[-1], t2[-1]), 1000)
y1 = linfit(xx, *params1)
y2 = linfit(xx, *params2)

plt.plot(xx, y1, color="red", label="Heizrate 1")
plt.plot(xx, y2, color="green", label="Heizrate 2")
print("Heizrate 1 über Fit = ", params1[0], " +- ", np.diag(np.sqrt(cov1))[0])
print("Heizrate 2 über Fit = ", params2[0], " +- ", np.diag(np.sqrt(cov2))[0])
plt.legend()
plt.grid("on")
plt.figure()

b1 = ufloat(params1[0], np.sqrt(np.diag(cov1))[0])
b2 = ufloat(params2[0], np.sqrt(np.diag(cov2))[0])

# Linearer Fit für den ersten Untergrund
fit_T1 = np.concatenate((T1[:10], T1[-16:]))  # Anfangs- und Endwerte für den Fit
fit_I1 = np.concatenate((I1[:10], I1[-16:]))

params_u1, cov_u1 = curve_fit(linfit, fit_T1, fit_I1)
xx_fit1 = np.linspace(fit_T1[0], fit_T1[-1], 1000)
yy_fit1 = linfit(xx_fit1, *params_u1)

print("Parameter für Untergrund = ", params_u1, " +- ", np.sqrt(np.diag(cov_u1)))

plt.errorbar(T1, I1, yerr=I1_err, color='orange', marker='x', markersize=9, linestyle="none", label='1. Messreihe')
plt.plot(T1, I1, color="orange", linewidth=3, alpha=0.3)
plt.plot(xx_fit1, yy_fit1, color='red', linestyle='-', label='Untergrund') 
plt.xlabel('T in K')
plt.ylabel('I in A')
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

plt.errorbar(T2, I2, yerr=I2_err, color='blue', marker='+', markersize=9, linestyle="none", label='2. Messreihe')
plt.plot(T2, I2, color="blue", linewidth=3, alpha=0.3)
plt.plot(xx_fit2, yy_fit2, color='green', linestyle='-', label='Untergrund ')
plt.xlabel('T in K')
plt.ylabel('I in A')
plt.title('Der Depolarisationsstrom für die zweite Heizrate')
plt.grid(True)
plt.legend()
plt.figure()


# Daten bereinigen
I1 = I1 - linfit(T1, *params_u1)
I2 = I2 - linfit(T2, *params_u2)

plt.errorbar(T1, I1, yerr=I1_err, color="orange", marker="x", markersize=9, linestyle="none", label="1. Messreihe, bereinigt")
plt.plot(T1, I1, color="orange", linewidth=3, alpha=0.3)
plt.errorbar(T2, I2, yerr=I2_err, color="blue", marker="+", markersize=9, linestyle="none", label="2. Messreihe, bereinigt")
plt.plot(T2, I2, color="blue", linewidth=3, alpha=0.3)
plt.xlabel("T in K")
plt.ylabel("I in A")
plt.grid("on")
plt.legend()
plt.title("Bereinigte Messreihen")
plt.figure()



# Für einen sinnvollen Fit werden nur die gaußähnlichen Daten zwischen den positiven Minima gewählt

T1 = T1[I1 > 0]
I1 = I1[I1 > 0]
I1 = I1[T1 > 240]
T1 = T1[T1 > 240]
I1 = I1[T1 < 280]
T1 = T1[T1 < 280]
I1_err = I1_err[0:len(I1)]

T2 = T2[I2 > 0]
I2 = I2[I2 > 0]
I2 = I2[T2 > 240]
T2 = T2[T2 > 240]
I2 = I2[T2 < 280]
T2 = T2[T2 < 280]
I2_err = I2_err[0:len(I2)]


plt.errorbar(T1, I1, yerr=I1_err, color="orange", marker="x", markersize=9, linestyle="none", label="1. Messreihe, bereinigt")
plt.plot(T1, I1, color="orange", linewidth=3, alpha=0.3)
plt.errorbar(T2, I2, yerr=I2_err, color="blue", marker="+", markersize=9, linestyle="none", label="2. Messreihe, bereinigt")
plt.plot(T2, I2, color="blue", linewidth=3, alpha=0.3)
plt.xlabel("T in K")
plt.ylabel("I in A")
plt.grid("on")
plt.legend()
plt.title("Bereinigte Messreihen, interessante Messpunkte")
plt.figure()

# Polarisationsansatz
kB = constants.Boltzmann

# Daten umformen für Fit
ln_I1 = np.log(I1)
inv_T1 = 1/(kB*T1)
ln_I2 = np.log(I2)
inv_T2 = 1/(kB*T2)

# Daten für Fit manuell auswählen
fit_I1 =  ln_I1[2:10]
fit_T1 = inv_T1[2:10]

fit_I2 = ln_I2[1:8]
fit_T2 = inv_T2[1:8]

params_W1, cov_W1 = curve_fit(linfit, fit_T1, fit_I1)
xx1 = np.linspace(fit_T1[0]+2e18, fit_T1[-1]-2e18)
yy_W1 = linfit(xx1, *params_W1)

params_W2, cov_W2 = curve_fit(linfit, fit_T2, fit_I2)
xx2 = np.linspace(fit_T2[0]+2e18, fit_T2[-1]-2e18)
yy_W2 = linfit(xx2, *params_W2)

print("Fit 1 = ", params_W1, " +- ", np.sqrt(np.diag(cov_W1)))
print("Fit 2 = ", params_W2, " +- ", np.sqrt(np.diag(cov_W2)))


plt.scatter(inv_T1, ln_I1, color="orange", marker="x", label="1. Messreihe, bereinigt")
plt.plot(inv_T1, ln_I1, color="orange", linewidth=3, alpha=0.3)
plt.scatter(fit_T1, fit_I1, color="brown", marker="x", label="Daten für Fit 1")
plt.scatter(inv_T2, ln_I2, color="blue", marker="+", label="2. Messreihe, bereinigt")
plt.plot(inv_T2, ln_I2, color="blue", linewidth=3, alpha=0.3)
plt.scatter(fit_T2, fit_I2, color="lightblue", marker="+", label="Daten für Fit 2")

plt.plot(xx1, yy_W1, color="red", label="Fit 1")
plt.plot(xx2, yy_W2, color="green", label="Fit 2")

plt.grid("on")
plt.xlabel("1/(kBT) in 1/J")
plt.ylabel("ln(I/(1A))")
plt.legend()
plt.figure()



W1 = ufloat(-params_W1[0], np.sqrt(np.diag(cov_W1))[0])
W2 = ufloat(-params_W2[0], np.sqrt(np.diag(cov_W2))[0])
W = (W1+W2)/2

T1_max = T1[I1 == max(I1)]
T2_max = T2[I2 == max(I2)]

t0_1 = kB*T1_max**2/(b1*W1) * unp.exp(-W1/(kB*T1_max))
t0_2 = kB*T2_max**2/(b2*W2) * unp.exp(-W2/(kB*T2_max))

print("Relaxationszeit 1 = ", t0_1)
print("Relaxationszeit 2 = ", t0_2)

W_1 = W1
t1_1 = t0_1
W_2 = W2
t1_2 = t0_2

##### Stromdichte
# Das Integral wird endlich bei t

t1, T1x, I1x = np.genfromtxt("Daten1.txt", unpack=True)
t2, T2x, I2x = np.genfromtxt("Daten2.txt", unpack=True)
I1x = I1x * 10**(-11)
T1x = T1x+273.15
I2x = I2x * 10**(-11)
T2x = T2x+273.15

# Zeiten des Gaußpeaks
t1 = t1[np.isin(T1x, T1)]
t2 = t2[np.isin(T2x, T2)]


# Idee fürs Integral = einfach Zeitschritt = 1 min 
def F(T, I):
    F_arr = np.zeros(len(T))

    for index in range(0, len(T)):
        integral = np.sum(I[index:-1])
        F = np.log(integral) - np.log(I[index])
        F_arr[index] = F
    return F_arr

F1 = F(T1, I1)
F2 = F(T2, I2)

T1 = T1[F1 >= 0]
F1 = F1[F1 >= 0]
T2 = T2[F2 >= 0]
F2 = F2[F2 >= 0]

T1_inv = 1/(kB * T1)
T2_inv = 1/(kB * T2)


params1, cov1 = curve_fit(linfit, T1_inv, F1)
params2, cov2 = curve_fit(linfit, T2_inv, F2)

xx = np.linspace(min(T1_inv[0], T2_inv[0])+2e18, max(T1_inv[-1], T2_inv[-1])-2e18)
yy1 = linfit(xx, *params1)
yy2 = linfit(xx, *params2)

W1 = ufloat(params1[0], np.sqrt(np.diag(cov1))[0])
t0_1 = unp.exp(ufloat(params1[1], np.sqrt(np.diag(cov1))[1]))

W2 = ufloat(params2[0], np.sqrt(np.diag(cov2))[0])
t0_2 = unp.exp(ufloat(params2[1], np.sqrt(np.diag(cov2))[1]))

print("W1 = ", W1)
print("t0_1 = ", t0_1)
print("W2 = ", W2)
print("t0_2 = ", t0_2)

W_3 = W1
t1_3 = t0_1
W_4 = W2
t1_4 = t0_2

plt.scatter(T1_inv, F1, color="orange", marker="x", label="1. Messreihe")
plt.plot(T1_inv, F1, color="orange", linewidth=3, alpha=0.3)
plt.scatter(T2_inv, F2, color="blue", marker="+", label="2. Messreihe")
plt.plot(T2_inv, F2, color="blue", linewidth=3, alpha=0.3)
plt.plot(xx, yy1, color="red", label="Fit 1")
plt.plot(xx, yy2, color="green", label="Fit 2")
plt.xlabel("1/kBT in 1/J")
plt.ylabel("F(T)")
plt.grid("on")
plt.legend()
plt.figure()



def relaxationszeit(T, W, t0):
    return unp.nominal_values(t0) * np.exp(unp.nominal_values(W)/(kB*T))

xx = np.linspace(270, 320)

t1 = relaxationszeit(xx, W_1, t1_1)
t2 = relaxationszeit(xx, W_2, t1_2)

plt.plot(xx, t1, color="orange", label="t, 1. Messreihe")
plt.plot(xx, t2, color="blue", label="t, 2. Messreihe")
plt.yscale("log")
plt.title("Relaxationszeit aus dem Polarisationsansatz")
plt.xlabel("T in K")
plt.ylabel("t in s")
plt.grid("on")
plt.figure()

t3 = relaxationszeit(xx, W_3, t1_3)
t4 = relaxationszeit(xx, W_4, t1_4)

plt.plot(xx, t3, color="orange", label="t, 1. Messreihe")
plt.plot(xx, t4, color="blue", label="t, 2. Messreihe")
plt.yscale("log")
plt.title("Relaxationszeit aus dem Stromdichtenansatz")
plt.xlabel("T in K")
plt.ylabel("t in s")
plt.grid("on")

plt.show()

"""
Heizrate 1 über Fit =  1.2420890575344499  +-  0.019833746560421136
Heizrate 2 über Fit =  1.6799985616114808  +-  0.019616136972691517

Parameter für Untergrund =  [ 0.0122142  -2.31501693]  +-  [0.00037518 0.09826862]
Parameter für Untergrund =  [ 0.02114334 -4.35632341]  +-  [0.00048343 0.12130725]

Polarisationsansatz
            -W                 C                    Delta W         Delta C
Fit 1 =  [-1.61487633e-19  1.92516720e+01]  +-  [1.18252214e-20 3.40514306e+00]
Fit 2 =  [-1.60155768e-19  1.88307760e+01]  +-  [1.20823475e-20 3.46349532e+00]

Relaxationszeit 1 =  [8.634802285151585e-20+/-2.934815305108716e-19]
Relaxationszeit 2 =  [1.4454749838472648e-19+/-4.975264060887047e-19]

Stromdichte
W1 =  (1.97+/-0.14)e-19
t0_1 =  (0.7+/-2.9)e-23
W2 =  (1.85+/-0.15)e-19
t0_2 =  (2+/-8)e-22


"""