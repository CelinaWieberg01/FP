import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# === Eingangsspannung ===
U_e = 0.107  # in Volt

# === Daten einlesen ===
f_all, U_a_all = np.genfromtxt("data/Invamp1.txt", unpack=True)

# === Verstärkung berechnen ===
V_all = U_a_all / U_e

# === Logarithmische Werte für alle Daten (nur für Plotzwecke) ===
ln_f_all = np.log(f_all)
ln_V_all = np.log(V_all)

# === Fit nur ab dem 9. Messwert durchführen ===
f_fit = f_all[8:]
V_fit = V_all[8:]

# === Potenzgesetz-Funktion ===
def power_law(f, a, b):
    return a * f**b

# === Direkter Fit auf V(f) = a * f^b ===
params, cov = curve_fit(power_law, f_fit, V_fit)
a, b = params
a_err, b_err = np.sqrt(np.diag(cov))

# === Plateau-Verstärkung bestimmen (aus allen Daten) ===
U_konst = np.mean(V_all[f_all < 1000])
ln_U_konst = np.log(U_konst)

# === Grenzfrequenz berechnen ===
f_cutoff = (U_konst / a)**(1 / b)
ln_f_cutoff = np.log(f_cutoff)

# === Logarithmierte Fit-Kurve für Plot ===
ln_f_fit = np.log(f_fit)
ln_V_fit_model = np.log(power_law(f_fit, a, b))



# === Fit nur ab dem 9. Messwert durchführen ===
f_fit = f_all[14:]
V_fit = V_all[14:]

# === Direkter Fit auf V(f) = a * f^b ===
params, cov = curve_fit(power_law, f_fit, V_fit)
a, b = params
a_err, b_err = np.sqrt(np.diag(cov))


# === Logarithmierte Fit-Kurve für Plot ===
ln_f_fit = f_fit
ln_V_fit_model = power_law(f_fit, a, b)

# === Plot erstellen ===
plt.figure(figsize=(8, 5))
plt.scatter(f_all, V_all, label="Alle Messdaten", color="blue")
plt.plot(ln_f_fit, ln_V_fit_model, label="Fit ab Punkt 9", color="red")
plt.axhline(ln_U_konst, color="green", linestyle="--", label=r"$U_{konst}$")
plt.axvline(ln_f_cutoff, color="gray", linestyle=":", label=f"$f_{{Grenze}}$ ≈ {ln_f_cutoff:.2f}")

# === Achsen und Legende ===
plt.xscale("log")
plt.yscale("log")
plt.xlabel("F [Hz]")
plt.ylabel(r"$V = \frac{U_{aus}}{U_{ein}}$")
plt.legend()
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# === Fit-Ergebnisse ausgeben ===
print(f"Fit-Parameter: a1 = {a:.4f} ± {a_err:.4f}, b1 = {b:.4f} ± {b_err:.4f}")
print(f"Gemittelte Plateau-Verstärkung U_konst1 = {U_konst:.4f}")
print(f"Grenzfrequenz f_cutoff1 ≈ {f_cutoff:.1f} Hz")

####################################################################################################################
# === Daten einlesen ===
f_all2, U_a_all2 = np.genfromtxt("data/Invamp2.txt", unpack=True)

# === Verstärkung berechnen ===
V_all2 = U_a_all2 / U_e

# === Potenzgesetz-Funktion ===
def power_law2(f2, a2, b2):
    return a2 * f2**b2

# === Fit durchführen ===
f_fit2 = f_all2[5:-2]
V_fit2 = V_all2[5:-2]
params2, cov2 = curve_fit(power_law2, f_fit2, V_fit2)
a2, b2 = params2
a_err2, b_err2 = np.sqrt(np.diag(cov2))

# === Plateau-Verstärkung bestimmen (aus allen Daten) ===
U_konst2 = np.mean(V_all2[f_all2 < 1000])

# === Grenzfrequenz berechnen ===
f_cutoff2 = (U_konst2 / a2)**(1 / b2)

# === Fit-Kurve für Plot ===
V_fit_model2 = power_law2(f_fit2, a2, b2)

# === Plot erstellen ===
plt.figure(figsize=(8, 5))
plt.scatter(f_all2, V_all2, label="Alle Messdaten", color="blue")
plt.plot(f_fit2, V_fit_model2, label="Fit", color="red")
plt.axhline(U_konst2, color="green", linestyle="--", label=r"$U_{konst}$")
plt.axvline(f_cutoff2, color="gray", linestyle=":", label=f"$f_{{Grenze}}$ ≈ {f_cutoff2:.1f} Hz")

# === Achsen und Legende ===
plt.xscale("log")
plt.yscale("log")
plt.xlabel("F [Hz]")
plt.ylabel(r"$V = \frac{U_{aus}}{U_{ein}}$")
plt.legend()
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
#plt.title("Messreihe 2: Fit ab Punkt 6")
plt.show()

# === Fit-Ergebnisse ausgeben ===
print(f"Fit-Parameter: a2 = {a2:.4f} ± {a_err2:.4f}, b2 = {b2:.4f} ± {b_err2:.4f}")
print(f"Gemittelte Plateau-Verstärkung U_konst2 = {U_konst2:.4f}")
print(f"Grenzfrequenz f_cutoff2 ≈ {f_cutoff2:.1f} Hz")

###############################################################################################
# === Daten einlesen ===
f_all3, U_a_all3 = np.genfromtxt("data/Invamp3.txt", unpack=True)

# === Verstärkung berechnen ===
V_all3 = U_a_all3 / U_e

# === Potenzgesetz-Funktion ===
def power_law3(f3, a3, b3):
    return a3 * f3**b3

# === Fit durchführen ===
f_fit3 = f_all3[8:]
V_fit3 = V_all3[8:]
params3, cov3 = curve_fit(power_law3, f_fit3, V_fit3)
a3, b3 = params3
a_err3, b_err3 = np.sqrt(np.diag(cov3))

# === Plateau-Verstärkung bestimmen (aus allen Daten) ===
U_konst3 = np.mean(V_all3[f_all3 < 1000])

# === Grenzfrequenz berechnen ===
f_cutoff3 = (U_konst3 / a3)**(1 / b3)

# === Fit-Kurve für Plot ===
V_fit_model3 = power_law3(f_fit3, a3, b3)

# === Plot erstellen ===
plt.figure(figsize=(8, 5))
plt.scatter(f_all3, V_all3, label="Alle Messdaten", color="blue")
plt.plot(f_fit3, V_fit_model3, label="Fit", color="red")
plt.axhline(U_konst3, color="green", linestyle="--", label=r"$U_{konst}$")
plt.axvline(f_cutoff3, color="gray", linestyle=":", label=f"$f_{{Grenze}}$ ≈ {f_cutoff3:.1f} Hz")

# === Achsen und Legende ===
plt.xscale("log")
plt.yscale("log")
plt.xlabel("F [Hz]")
plt.ylabel(r"$V = \frac{U_{aus}}{U_{ein}}$")
plt.legend()
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
#plt.title("Messreihe 2: Fit ab Punkt 6")
plt.show()

# === Fit-Ergebnisse ausgeben ===
print(f"Fit-Parameter: a3 = {a3:.4f} ± {a_err3:.4f}, b3 = {b3:.4f} ± {b_err3:.4f}")
print(f"Gemittelte Plateau-Verstärkung U_konst3 = {U_konst3:.4f}")
print(f"Grenzfrequenz f_cutoff3 ≈ {f_cutoff3:.1f} Hz")
############################################################################
