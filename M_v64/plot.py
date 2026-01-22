import matplotlib
#matplotlib.use("MacOSX")

#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import sem

# Dateien laden

file1 = np.genfromtxt("data/contrast1.txt", comments="#")
file2 = np.genfromtxt("data/contrast2.txt", comments="#")
file3 = np.genfromtxt("data/contrast3.txt", comments="#")

angles = file1[:, 0]
Imax_all = np.vstack([file1[:, 1], file2[:, 1], file3[:, 1]])
Imin_all = np.vstack([file1[:, 2], file2[:, 2], file3[:, 2]])

#  Mittelwerte und Standardabweichungen

Imax_mean = np.mean(Imax_all, axis=0)
Imax_std = sem(Imax_all, axis=0)

Imin_mean = np.mean(Imin_all, axis=0)
Imin_std = sem(Imin_all, axis=0)



#  Kontrast berechnen 

Imax_u = unp.uarray(Imax_mean, Imax_std)
Imin_u = unp.uarray(Imin_mean, Imin_std)

V_u = (Imax_u - Imin_u) / (Imax_u + Imin_u)

# Betrag der Messdaten
V = np.abs(unp.nominal_values(V_u))
V_err = unp.std_devs(V_u)


# Fit-Modell

def contrast_model(phi, Voff, V0):
    return Voff + V0 * np.abs(np.sin(2 * np.deg2rad(phi)))

popt, pcov = curve_fit(
    contrast_model,
    angles,
    V,
    sigma=V_err,
    absolute_sigma=True
)

Voff_fit, V0_fit = popt
Voff_err, V0_err = np.sqrt(np.diag(pcov))

print("Fit-Ergebnisse:")
print(f"V0 = {V0_fit:.4f} ± {V0_err:.4f}")
print(f"Voff = {Voff_fit:.4f} ± {Voff_err:.4f}")


#  Plot

phi_plot = np.linspace(0, 180, 500)

plt.figure(figsize=(8, 5))

# Messdaten mit Fehlerbalken
plt.errorbar(
    angles, V, yerr=V_err,
    fmt='o', markersize=6, capsize=4,
    label=r"calculated contrast $C$", color="black"
)

# Fitkurve
plt.plot(
    phi_plot,
    contrast_model(phi_plot, *popt),
    color="green",
    linewidth=2,
    label=r"Fit"
)

plt.xlabel(r"Polarization angle $\phi$ in $\si{\degree}$", fontsize=14)

plt.ylabel(r"$C(\phi)$", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.title(r"Contrast at different polarization angles")
plt.savefig("plots/contrast.pdf")
plt.figure()

########################################################## Glas ########################################################

# Gegebene Konstanten 

lambda_vac = 632.99e-9 
d = 1e-3 
theta0_deg = 10 
theta0 = np.deg2rad(theta0_deg) 

# Unsicherheiten
dM = 2 # Zählunsicherheit
dtheta_deg = 2 # Winkelunsicherheit 
dtheta = np.deg2rad(dtheta_deg)

theta = ufloat(theta0, dtheta)
# Tabelle der gemessenen Maxima

M_values = np.array([29, 36, 32, 25, 32, 34, 27, 32, 32, 33, 32, 29])
M_errs = 2*np.ones(len(M_values))
M_u = unp.uarray(M_values, M_errs)

# Formel für n(M)
def refractive_index(M, theta): return 1 / (1 - (lambda_vac * M) / (2 * d * theta0 * theta))



# n für jede Messung 

n_values = refractive_index(M_values, theta)


# Mittelwert und Standardabweichung

n_mean = np.mean(n_values) 


print("Einzelne n-Werte:")
print(n_values) 

print("\nMittelwert des Brechungsindex:")
print(f"n_mean = {n_mean:.4f}") 

######################################################################################################################
# 
# # 1. Konstanten
# 
# lambda_vac = 632.996e-9        
# L = 100.0e-3                   
# dL = 0.1e-3                    
# 
# R = 8.3144                     
# T = 22.0 + 273.15              
# 
# # Zielwerte
# T0 = 15.0 + 273.15             
# p0 = 1013.0 * 100              
# 
# # 2. Drei Messreihen laden
# 
# file1 = np.genfromtxt("data/air1.txt")
# file2 = np.genfromtxt("data/air2.txt")
# file3 = np.genfromtxt("data/air3.txt")
# 
# # Jede Datei: p[mbar], counts
# p1, M1 = file1[:,0], file1[:,1]
# p2, M2 = file2[:,0], file2[:,1]
# p3, M3 = file3[:,0], file3[:,1]
# 
# 
# # 3. Mittelwerte bilden
# 
# p_mbar = p1   #
# p = p_mbar * 100   
# 
# M_mean = (M1 + M2 + M3) / 3.0
# 
# 
# # 4. Δn(p) berechnen
# 
# delta_n = M_mean * lambda_vac / L
# 
# # 5. Linearer Fit Δn(p) 
# 
# def linear(p, c, b):
#     return c*p + b
# 
# popt, pcov = curve_fit(linear, p, delta_n)
# c_fit, b_fit = popt
# c_err, b_err = np.sqrt(np.diag(pcov))
# 
# print("\nFit-Ergebnisse:")
# print(f"c = {c_fit:.3e} ± {c_err:.3e}   Pa^-1")
# print(f"b = {b_fit:.3e} ± {b_err:.3e}")
# 
# 
# # 6. Molar refractivity 
# 
# A = c_fit * (2 * R * T) / 3
# A_err = c_err * (2 * R * T) / 3
# 
# print("\nMolar refractivity A:")
# print(f"A = {A:.3e} ± {A_err:.3e}   m^3/mol")
# 
# # 7. Refractive index 
# 
# delta_n_T0_p0 = (A * 3 / (2 * R * T0)) * p0 + b_fit
# delta_n_T0_p0_err = np.sqrt(
#     (A_err * 3/(2*R*T0) * p0)**2 +
#     b_err**2
# )
# 
# n_air = 1 + delta_n_T0_p0
# n_air_err = delta_n_T0_p0_err
# 
# print("\nRefractive index at T0, p0:")
# print(f"Δn_air = {delta_n_T0_p0:.7f} ± {delta_n_T0_p0_err:.7f}")
# print(f"n_air  = {n_air:.7f} ± {n_air_err:.7f}")
# 
# 
# # 8. Plot
# 
# p_plot = np.linspace(min(p), max(p), 500)
# 
# plt.figure(figsize=(8,5))
# plt.plot(p_plot, linear(p_plot, *popt), 'g-', label="Fit")
# plt.scatter(p, delta_n, color='black', label="Data")
# 
# plt.xlabel("p [pa]")
# plt.ylabel(r"$\Delta n$")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()


Lambda_vac = 632.996e-9 
L = 100.0e-3                   
dL = 0.1e-3
L_u = ufloat(L, dL)

def linfit(x, m, b):
    return m*x+b

def delta_n_func(M):
    return M*Lambda_vac/L_u 

p1, M1 = np.genfromtxt("data/air1.txt", unpack=True)
p2, M2 = np.genfromtxt("data/air2.txt", unpack=True)
p3, M3 = np.genfromtxt("data/air3.txt", unpack=True)

p1 *= 100
p2 *= 100
p3 *= 100

delta_n1 = delta_n_func(M1)
delta_n2 = delta_n_func(M2)
delta_n3 = delta_n_func(M3)

p = np.concatenate((p1, p2, p3))
delta_n = np.concatenate((delta_n1, delta_n2, delta_n3))

params, cov = curve_fit(linfit, p, unp.nominal_values(delta_n))
m = ufloat(params[0], np.sqrt(np.diag(cov))[0])
b = ufloat(params[1], np.sqrt(np.diag(cov))[1])

print("fit m*x+b für delta n in p in mbar:")
print(f"m = {m:.4}")
print(f"b = {b:.4}")

T = 22.0 + 273 # K
R = 8.3144 # J/molK
A = m*2*R*T/3
print(f"A = {A}")

n_standard = A * 3/2/R * (101300)/(15+273) + b + 1
print(n_standard)

plt.errorbar(p1, unp.nominal_values(delta_n1), yerr=unp.std_devs(delta_n1), fmt=".", label=r"Measurement 1")
plt.errorbar(p2, unp.nominal_values(delta_n2), yerr=unp.std_devs(delta_n2), fmt=".", label=r"Measurement 2")
plt.errorbar(p3, unp.nominal_values(delta_n3), yerr=unp.std_devs(delta_n3), fmt=".", label=r"Measurement 3")

xx = np.linspace(0, 102000)
plt.plot(xx, linfit(xx, *params), label="Fit")

plt.legend()
plt.xlabel(r"$p$ in $\si{\pascal}$")
plt.ylabel(r"$\Delta n$")
plt.title(r"Refractive index difference between vacuum and air")

plt.savefig("plots/air.pdf")