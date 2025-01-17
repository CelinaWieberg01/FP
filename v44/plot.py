import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl


#### DETEKROSCAN A)

alpha_detectorscan, intensity_detectorscan = np.genfromtxt("data/Detector1.UXD", unpack=True)

def gauss(x, my, sigma, a, b):
    return a/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((x-my)/sigma)**2) + b

params_detectorscan, cov_detectorscan = curve_fit(gauss, alpha_detectorscan, intensity_detectorscan, p0=[0, 3, 350000, 10000])

xx = np.linspace(-0.5, 0.5, 1000)
yy = gauss(xx, *params_detectorscan)

fullmaximum = params_detectorscan[2]/(params_detectorscan[1]*np.sqrt(2*np.pi))
halfmaximum = fullmaximum / 2
halfwidth = xx[yy >= halfmaximum]


plt.plot(xx, gauss(xx, *params_detectorscan), color="red", label="Fit")
plt.hlines(halfmaximum, halfwidth[0], halfwidth[-1], "black", label="HWFM")
plt.scatter(alpha_detectorscan, intensity_detectorscan, marker="x", s=10, label="Daten")
plt.grid("on")
plt.legend()
plt.xlabel("alpha")
plt.ylabel("Intensität")
plt.savefig("plots/detectorscan.pdf")
plt.figure()

# print("Detektorscan:")
# print("my = ", params_detectorscan[0], "+-", np.sqrt(np.diag(cov_detectorscan))[0])
# print("sigma = ", params_detectorscan[1], "+-", np.sqrt(np.diag(cov_detectorscan))[1])
# print("a = ", params_detectorscan[2], "+-", np.sqrt(np.diag(cov_detectorscan))[2])
# print("b = ", params_detectorscan[3], "+-", np.sqrt(np.diag(cov_detectorscan))[3])
# print("fullmaximum = ", fullmaximum)
# print("halfmaximum = ", halfmaximum)
# print("halfwidth = ", halfwidth[0], " - ", halfwidth[-1])

"""
Detektorscan
my =  0.002676671441559556 +- 0.00046026460585211495
sigma =  0.03696928070868331 +- 0.0004782276504064794
a =  34247.137989212286 +- 410.97008376851716
b =  5197.864799863819 +- 917.7848278672068

fullmaximum =  369567.13965575804
halfmaximum =  184783.56982787902
halfwidth =  -0.04154154154154155  -  0.04654654654654655
"""








##############################################

### GEOMETRIEFAKTOR C)


# ERST STRAHLBREITE

z, intensity_z = np.genfromtxt("data/Zscan1.UXD", unpack=True)

mittlere_full_intensity_z = np.mean(intensity_z[z <= 0])
mittlere_no_intensity_z = np.mean(intensity_z[intensity_z <= 0.01 * mittlere_full_intensity_z])

z_wafer = z[np.logical_and(intensity_z > mittlere_no_intensity_z, intensity_z < mittlere_full_intensity_z)]

intensity_wafer = intensity_z[np.logical_and(intensity_z > mittlere_no_intensity_z, intensity_z < mittlere_full_intensity_z)]
intensity_wafer = intensity_wafer[np.logical_and(z_wafer > 0.25, z_wafer < 0.5)]

z_wafer = z_wafer[np.logical_and(z_wafer > 0.25, z_wafer < 0.5)]

def linfit(x, m, b):
    return m*x + b

### FIT FÜR WAFER
params_wafer, cov_wafer = curve_fit(linfit, z_wafer, intensity_wafer, p0=[-14000000, 700000])
xx = np.linspace(0.29, 0.48, 10000)
yy = linfit(xx, *params_wafer)

strahlhöhe_unten = xx[yy < mittlere_full_intensity_z][0]
strahlhöhe_oben = xx[yy > mittlere_no_intensity_z][-1]

strahlbreite = strahlhöhe_oben - strahlhöhe_unten

plt.hlines(mittlere_full_intensity_z, -1, 0.4, colors="red", linestyles="dashed", label="Maixmale Intensität")
plt.hlines(mittlere_no_intensity_z, 0.4, 1, colors="green", linestyles="dashed", label="Minimale Intensität")
plt.vlines(strahlhöhe_unten, 0, mittlere_full_intensity_z, color="grey", linestyles="dashed", label="Strahlgrenzen")
plt.vlines(strahlhöhe_oben, 0, mittlere_full_intensity_z, color="grey", linestyles="dashed")

plt.plot(xx, yy, color="lightblue", label="Fit")

plt.scatter(z_wafer, intensity_wafer, s=30, marker="o", label="Wafer",color="orange")
plt.scatter(z, intensity_z, s=10, marker="x", label="Daten", color="blue")

plt.xlabel("z in mm")
plt.ylabel("Intensität")
plt.grid("on")
plt.legend()
plt.savefig("plots/zscan.pdf")
plt.figure()

print("Strahlbreite = ", strahlbreite)
print("mittlere volle Intensität = ", mittlere_full_intensity_z)
print("mittlere no intensität = ", mittlere_no_intensity_z)
print("Fitparameter:")
print("a = ", params_wafer[0], " +- ", np.sqrt(np.diag(cov_wafer))[0])
print("b = ", params_wafer[1], " +- ", np.sqrt(np.diag(cov_wafer))[1])
print("Lasergrenzen laut Fit = ", strahlhöhe_unten, " und ", strahlhöhe_oben, ". Distanz = ", strahlbreite)




alpha_rocking, intensity_rocking = np.genfromtxt("data/Rocking1_2.UXD", unpack=True)
alpha_rocking_rechts = alpha_rocking[alpha_rocking > 0]
intensity_rocking_rechts = intensity_rocking[alpha_rocking > 0]
geometriewinkel_rechts = alpha_rocking_rechts[intensity_rocking_rechts <= 0.01 * max(intensity_rocking)][0]
alpha_rocking_links = alpha_rocking[alpha_rocking < 0]
intensity_rocking_links = intensity_rocking[alpha_rocking < 0]
geometriewinkel_links = alpha_rocking_links[intensity_rocking_links <= 0.01 * max(intensity_rocking)][-1]
mittlerer_geometriewinkel = (abs(geometriewinkel_links) + abs(geometriewinkel_rechts))/2
abweichung_mittlerer_geometriewinkel = 0.5 * ( ( abs(abs(geometriewinkel_links) - mittlerer_geometriewinkel) + abs(abs(geometriewinkel_rechts) - mittlerer_geometriewinkel)))

D = 20 # mm, länge des wafers

def radtodeg(rad):
    return rad*57.2958

def degtorad(deg):
    return deg*np.pi/180

alpha_g_theo = unp.arcsin(ufloat(strahlbreite, 0.030)/D)
alpha_g_theo_nom = radtodeg(unp.nominal_values(alpha_g_theo))
alpha_g_theo_std = radtodeg(unp.std_devs(alpha_g_theo))
print("Geometriewinkel-Theorie = ", ufloat(alpha_g_theo_nom, alpha_g_theo_std))

plt.vlines([geometriewinkel_links, geometriewinkel_rechts], 0, max(intensity_rocking), linestyles="dashed", color="red")

plt.scatter(alpha_rocking, intensity_rocking, s=10, marker="x", label="Daten")
plt.grid("on")
plt.legend()
plt.xlabel("alpha")
plt.ylabel("intensität")
plt.savefig("plots/rockingscan.pdf")
plt.figure()

print("Geometriewinkel = ", geometriewinkel_links, " und ", geometriewinkel_rechts)
print("mittlerer Geometriewinkel = ", mittlerer_geometriewinkel, " +- ", abweichung_mittlerer_geometriewinkel)


"""
Strahlbreite =  0.17631863186318636
mittlere volle Intensität =  374912.26923076925
mittlere no intensität =  399.6666666666667
Fitparameter:
a =  -2123867.4963938277  +-  241153.8795209772
b =  1002844.098701778  +-  87880.65573862198
Lasergrenzen laut Fit =  0.2956625662566256  und  0.471981198119812 . Distanz =  0.17631863186318636
Geometriewinkel-Theorie =  0.505122396592914
Geometriewinkel =  -0.4  und  0.44
mittlerer Geometriewinkel =  0.42000000000000004  +-  0.01999999999999999
"""




#######################################

# REFLEKTIVITÄTSSCAN DIFFUSE SCAN B)

alpha_reflect, intensity_reflect = np.genfromtxt("data/Reflect1.UXD", unpack=True)
alpha_diffuse, intensity_diffuse = np.genfromtxt("data/Reflect2.UXD", unpack=True)

mittlerer_geometriewinkel_rad = degtorad(mittlerer_geometriewinkel)

def geometrische_korrektur_reflexion(winkel, intensity):
    winkel = degtorad(winkel)
    G = (D/strahlbreite * np.sin(winkel[winkel < mittlerer_geometriewinkel_rad]))
    intensity[winkel < mittlerer_geometriewinkel_rad] /= G
    print(G)
    return intensity

smaller_length = min(len(alpha_diffuse), len(alpha_reflect))
intensity_true = intensity_reflect[:smaller_length] - intensity_diffuse[:smaller_length]
alpha_true = alpha_reflect[:smaller_length]

alpha_geometrische_korrektur = alpha_true
intensity_geometrische_korrektur = geometrische_korrektur_reflexion(alpha_geometrische_korrektur, intensity_true)

# Wellenlänge
wavelength = 1.541*10**(-10) # m
r_e = 2.817*10**(-15) # m 
rho = 1.08 * 10**10 / (10**(-6))
critical_angle = wavelength * np.sqrt(r_e * rho / np.pi)
# ??????

def ideal_fresnel(alpha_i):
    R_F = (0.223/(2*alpha_i + 1e-8))**4
    return R_F

testangles = np.linspace(0, alpha_reflect[-1], 1000)
ideal_fresnel_reflect = ideal_fresnel(testangles)

plt.plot(alpha_reflect, intensity_reflect, label="Reflektivitätsscan")
plt.plot(alpha_diffuse, intensity_diffuse, label="Diffuser Scan")
plt.plot(alpha_true,    intensity_true,    label="Differenz der Scans")
plt.plot(alpha_geometrische_korrektur, intensity_geometrische_korrektur, label="Geometrische Korrektur")
plt.grid("on")
plt.legend()
plt.yscale("log")
plt.xlabel("alpha")
plt.ylabel("Intensität")
plt.savefig("plots/reflectivity.pdf")
plt.figure()

