import matplotlib
matplotlib.use("MacOSX")

matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit


# Dateien laden

file1 = np.genfromtxt("data/contrast1.txt", comments="#")
file2 = np.genfromtxt("data/contrast2.txt", comments="#")
file3 = np.genfromtxt("data/contrast3.txt", comments="#")

angles = file1[:, 0]
Imax_all = np.vstack([file1[:, 1], file2[:, 1], file3[:, 1]])
Imin_all = np.vstack([file1[:, 2], file2[:, 2], file3[:, 2]])

#  Mittelwerte und Standardabweichungen

Imax_mean = np.mean(Imax_all, axis=0)
Imax_std = np.std(Imax_all, axis=0, ddof=1)

Imin_mean = np.mean(Imin_all, axis=0)
Imin_std = np.std(Imin_all, axis=0, ddof=1)



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
    label=r"Data ($|V(\phi)|$)", color="black"
)

# Fitkurve
plt.plot(
    phi_plot,
    contrast_model(phi_plot, *popt),
    color="green",
    linewidth=2,
    label=r"Fit"
)

plt.xlabel(r"Polarization angle $\phi$ ", fontsize=14)

plt.ylabel(r"$V(\phi)$", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)

plt.show()

########################################################## Glas ########################################################

# Gegebene Konstanten 

lambda_vac = 632.99e-9 
d = 1e-3 
theta0_deg = 10 
theta0 = np.deg2rad(theta0_deg) 

# Unsicherheiten
dM = 1 # Zählunsicherheit
dtheta_deg = 1 # Winkelunsicherheit 
dtheta = np.deg2rad(dtheta_deg)

# Tabelle der gemessenen Maxima

M_values = np.array([29, 36, 32, 25, 32, 34, 27, 32, 32, 33, 32, 29]) 

# Formel für n(M)
def refractive_index(M, theta): return 1 / (1 - (lambda_vac * M) / (2 * d * theta0 * theta))

theta = theta0

# n für jede Messung 

n_values = refractive_index(M_values, theta)

# Mittelwert und Standardabweichung

n_mean = np.mean(n_values) 
n_std = np.std(n_values, ddof=1) 


dn_dM = (lambda_vac / (2 * d * theta0 * theta)) / (1 - (lambda_vac * M_values) / (2 * d * theta0 * theta))**2 
dn_dtheta = (lambda_vac * M_values / (2 * d * theta0 * theta**2)) / (1 - (lambda_vac * M_values) / (2 * d * theta0 * theta))**2

# Gesamtfehler pro Messung
n_errors = np.sqrt((dn_dM * dM)**2 + (dn_dtheta * dtheta)**2)

# Gemittelter Fehler
n_error_mean = np.sqrt(np.sum(n_errors**2)) / len(n_errors)

print("Einzelne n-Werte:")
print(n_values) 

print("\nMittelwert des Brechungsindex:")
print(f"n_mean = {n_mean:.4f}") 

print("\nStandardabweichung:") 
print(f"n_std = {n_std:.4f}") 

print("\nFehler des Mittelwerts (inkl. Fehlerfortpflanzung):") 
print(f"n_error = {n_error_mean:.4f}")

print("\nErgebnis:") 
print(f"n_glass = {n_mean:.4f} ± {n_error_mean:.4f}")

######################################################################################################################

# 1. Konstanten

lambda_vac = 632.996e-9        
L = 100.0e-3                   
dL = 0.1e-3                    

R = 8.3144                     
T = 22.2 + 273.15              

# Zielwerte
T0 = 15.0 + 273.15             
p0 = 1013.0 * 100              

# 2. Drei Messreihen laden

file1 = np.genfromtxt("data/air1.txt")
file2 = np.genfromtxt("data/air2.txt")
file3 = np.genfromtxt("data/air3.txt")

# Jede Datei: p[mbar], counts
p1, M1 = file1[:,0], file1[:,1]
p2, M2 = file2[:,0], file2[:,1]
p3, M3 = file3[:,0], file3[:,1]


# 3. Mittelwerte bilden

p_mbar = p1   #
p = p_mbar * 100   

M_mean = (M1 + M2 + M3) / 3.0


# 4. Δn(p) berechnen

delta_n = M_mean * lambda_vac / L

# 5. Linearer Fit Δn(p) 

def linear(p, c, b):
    return c*p + b

popt, pcov = curve_fit(linear, p, delta_n)
c_fit, b_fit = popt
c_err, b_err = np.sqrt(np.diag(pcov))

print("\nFit-Ergebnisse:")
print(f"c = {c_fit:.3e} ± {c_err:.3e}   Pa^-1")
print(f"b = {b_fit:.3e} ± {b_err:.3e}")


# 6. Molar refractivity 

A = c_fit * (2 * R * T) / 3
A_err = c_err * (2 * R * T) / 3

print("\nMolar refractivity A:")
print(f"A = {A:.3e} ± {A_err:.3e}   m^3/mol")

# 7. Refractive index 

delta_n_T0_p0 = (A * 3 / (2 * R * T0)) * p0 + b_fit
delta_n_T0_p0_err = np.sqrt(
    (A_err * 3/(2*R*T0) * p0)**2 +
    b_err**2
)

n_air = 1 + delta_n_T0_p0
n_air_err = delta_n_T0_p0_err

print("\nRefractive index at T0, p0:")
print(f"Δn_air = {delta_n_T0_p0:.7f} ± {delta_n_T0_p0_err:.7f}")
print(f"n_air  = {n_air:.7f} ± {n_air_err:.7f}")


# 8. Plot

p_plot = np.linspace(min(p), max(p), 500)

plt.figure(figsize=(8,5))
plt.plot(p_plot, linear(p_plot, *popt), 'g-', label="Fit")
plt.scatter(p, delta_n, color='black', label="Data")

plt.xlabel("p [pa]")
plt.ylabel(r"$\Delta n$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
