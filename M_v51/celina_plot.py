import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# === Eingangsspannung ===
U_e = 2.5  # in Volt; passe für jede Messreihe an

# === Daten einlesen ===
f_all, U_a_all = np.genfromtxt("data/Invamp1.txt", unpack=True)

# === Verstärkung berechnen ===
V_all = U_a_all / U_e

# === Potenzgesetz-Funktion ===
def power_law(f, a, b):
    return a * f**b

# === Erster Fit (linker Abfall) ===
f_fit1 = f_all[1:8]      # Beispiel: setze hier sinnvolle Indizes
V_fit1 = V_all[1:8]
p1, cov1 = curve_fit(power_law, f_fit1, V_fit1, maxfev=10000)
a1, b1 = p1
a1_err, b1_err = np.sqrt(np.diag(cov1))

# === Zweiter Fit (rechter Abfall) ===
f_fit2 = f_all[8:]      # Beispiel: setze hier sinnvolle Indizes
V_fit2 = V_all[8:]
p2, cov2 = curve_fit(power_law, f_fit2, V_fit2, maxfev=10000)
a2, b2 = p2
a2_err, b2_err = np.sqrt(np.diag(cov2))

# === Plateau-Verstärkung ===
U_konst = np.mean(V_all[f_all < 1000])

# === Grenzfrequenzen ===
f_cut_left = (U_konst / a1)**(1.0 / b1)
f_cut_right = (U_konst / a2)**(1.0 / b2)
bandwidth = f_cut_right - f_cut_left

# === Plot (log-log) - benutze f und V, nicht ln-Variablen ===
plt.figure(figsize=(8,5))
plt.scatter(f_all, V_all, label="Messdaten")
f_plot = np.logspace(np.log10(f_all.min()*0.8), np.log10(f_all.max()*1.2), 300)
plt.plot(f_plot, power_law(f_plot, a1, b1), 'r-', label="Fit links")
plt.plot(f_plot, power_law(f_plot, a2, b2), 'b-', label="Fit rechts")
#plt.axhline(U_konst, color='green', linestyle='--', label='U_konst')
plt.axvline(f_cut_left, color='purple', linestyle=':', label=f'f_cut_links ≈ {f_cut_left:.1f} Hz')
plt.axvline(f_cut_right, color='purple', linestyle=':', label=f'f_cut_rechts ≈ {f_cut_right:.1f} Hz')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('F [Hz]'); plt.ylabel(r'V =$\frac{U_{aus}}{U_{ein}}$')
plt.legend(); plt.grid(True); plt.show()

print("a1,b1:", a1, b1, "a2,b2:", a2, b2)
print("U_konst, f_cut_left, f_cut_right, bandwidth:", U_konst, f_cut_left, f_cut_right, bandwidth)

####################################################################################################################

# === Daten einlesen ===
f_all2, U_a_all2 = np.genfromtxt("data/Invamp2.txt", unpack=True)

# === Verstärkung berechnen ===
V_all2 = U_a_all2 / U_e

# === Potenzgesetz-Funktion ===
def power_law2(f2, a2, b2):
    return a2 * f2**b2

# === Erster Fit (linker Abfall) ===
f_fit3 = f_all2[1:5]      # Beispiel: setze hier sinnvolle Indizes
V_fit3 = V_all2[1:5]
p3, cov3 = curve_fit(power_law2, f_fit3, V_fit3, maxfev=10000)
a3, b3 = p3
a3_err, b3_err = np.sqrt(np.diag(cov3))

# === Zweiter Fit (rechter Abfall) ===
f_fit4 = f_all2[7:-1]      # Beispiel: setze hier sinnvolle Indizes
V_fit4 = V_all2[7:-1]
p4, cov4 = curve_fit(power_law2, f_fit4, V_fit4, maxfev=10000)
a4, b4 = p4
a4_err, b4_err = np.sqrt(np.diag(cov4))

# === Plateau-Verstärkung ===
U_konst = np.mean(V_all2[f_all2 < 1000])

# === Grenzfrequenzen ===
f_cut_left2 = (U_konst / a3)**(1.0 / b3)
f_cut_right2 = (U_konst / a4)**(1.0 / b4)
bandwidth2 = f_cut_right2 - f_cut_left2

# === Plot (log-log) - benutze f und V, nicht ln-Variablen ===
plt.figure(figsize=(8,5))
plt.scatter(f_all2, V_all2, label="Messdaten")
f_plot = np.logspace(np.log10(f_all.min()*0.8), np.log10(f_all.max()*1.2), 300)
plt.plot(f_plot, power_law2(f_plot, a3, b3), 'r-', label="Fit links")
plt.plot(f_plot, power_law2(f_plot, a4, b4), 'b-', label="Fit rechts")
#plt.axhline(U_konst, color='green', linestyle='--', label='U_konst')
plt.axvline(f_cut_left2, color='purple', linestyle=':', label=f'f_cut_links ≈ {f_cut_left2:.1f} Hz')
plt.axvline(f_cut_right2, color='purple', linestyle=':', label=f'f_cut_rechts ≈ {f_cut_right2:.1f} Hz')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('F [Hz]'); plt.ylabel(r'V =$\frac{U_{aus}}{U_{ein}}$')
plt.legend(); plt.grid(True); plt.show()

print("a3,b3:", a3, b3, "a4,b4:", a4, b4)
print("U_konst, f_cut_left2, f_cut_right2, bandwidth2:", U_konst, f_cut_left2, f_cut_right2, bandwidth2)

####################################################################################################################
# === Daten einlesen ===
f_all3, U_a_all3 = np.genfromtxt("data/Invamp3.txt", unpack=True)

# === Verstärkung berechnen ===
V_all3 = U_a_all3 / U_e

# === Potenzgesetz-Funktion ===
def power_law3(f3, a3, b3):
    return a3 * f3**b3

# === Erster Fit (linker Abfall) ===
f_fit5 = f_all3[1:8]      # Beispiel: setze hier sinnvolle Indizes
V_fit5 = V_all3[1:8]
p5, cov5 = curve_fit(power_law3, f_fit5, V_fit5, maxfev=10000)
a5, b5 = p5
a5_err, b5_err = np.sqrt(np.diag(cov5))

# === Zweiter Fit (rechter Abfall) ===
f_fit6 = f_all3[8:]      # Beispiel: setze hier sinnvolle Indizes
V_fit6 = V_all3[8:]
p6, cov6 = curve_fit(power_law3, f_fit6, V_fit6, maxfev=10000)
a6, b6 = p6
a6_err, b6_err = np.sqrt(np.diag(cov6))

# === Plateau-Verstärkung ===
U_konst = np.mean(V_all3[f_all3 < 1000])

# === Grenzfrequenzen ===
f_cut_left3 = (U_konst / a5)**(1.0 / b5)
f_cut_right3 = (U_konst / a6)**(1.0 / b6)
bandwidth3 = f_cut_right3 - f_cut_left3

# === Plot (log-log) - benutze f und V, nicht ln-Variablen ===
plt.figure(figsize=(8,5))
plt.scatter(f_all3, V_all3, label="Messdaten")
f_plot = np.logspace(np.log10(f_all.min()*0.8), np.log10(f_all.max()*1.2), 300)
plt.plot(f_plot, power_law3(f_plot, a5, b5), 'r-', label="Fit links")
plt.plot(f_plot, power_law3(f_plot, a6, b6), 'b-', label="Fit rechts")
#plt.axhline(U_konst, color='green', linestyle='--', label='U_konst')
plt.axvline(f_cut_left3, color='purple', linestyle=':', label=f'f_cut_links ≈ {f_cut_left3:.1f} Hz')
plt.axvline(f_cut_right3, color='purple', linestyle=':', label=f'f_cut_rechts ≈ {f_cut_right3:.1f} Hz')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('F [Hz]'); plt.ylabel(r'V =$\frac{U_{aus}}{U_{ein}}$')
plt.legend(); plt.grid(True); plt.show()

print("a5,b5:", a5, b5, "a6,b6:", a6, b6)
print("U_konst, f_cut_left3, f_cut_right3, bandwidth3:", U_konst, f_cut_left3, f_cut_right3, bandwidth3)
################################################################################################################

def power_law7(f, a, b):
    return a * f**b

# -----------------------------
# === Integrator ===
# -----------------------------
# === Eingangsspannung  ===
U_e_int = 0.103  # Volt

# Parameter: C = 100 nF, R = 10 kΩ
C_int = 100e-9  # 100 nF
R_int = 10e3    # 10 kΩ

# === Daten einlesen ===
f_int, U_a_int = np.genfromtxt("data/IntegreatorAmp.txt", unpack=True)

# === Verstärkung berechnen ===
V_int = U_a_int / U_e_int

# === Nonlinearer Fit mit curve_fit ===
p_int, cov_int = curve_fit(power_law7, f_int, V_int, maxfev=10000)
a_int, b_int = p_int
a_err_int, b_err_int = np.sqrt(np.diag(cov_int))

# === Lineare Ausgleichsgerade im ln-ln Raum ===
ln_f_int = np.log(f_int)
ln_V_int = np.log(V_int)
coef_lin_int = np.polyfit(ln_f_int, ln_V_int, deg=1)  # coef_lin_int[0] = b_lin, coef_lin_int[1] = ln(a_lin)
b_lin_int = coef_lin_int[0]
a_lin_int = np.exp(coef_lin_int[1])

# === Zeitkonstante berechnen (wie in deinem Ansatz) ===
Ckonst_int = a_int * R_int * C_int

# === Plot Integrator: Daten, curve_fit (rot), ln-ln Ausgleich (orange gestrichelt) ===
plt.figure(figsize=(8,5))
plt.scatter(f_int, V_int, label="Integrator Messdaten", color="blue", s=20)

f_plot_int = np.logspace(np.log10(f_int.min()*0.8), np.log10(f_int.max()*1.2), 400)
plt.plot(f_plot_int, power_law7(f_plot_int, a_int, b_int), color="red", linewidth=2, label="Fit")
plt.plot(f_plot_int, power_law7(f_plot_int, a_lin_int, b_lin_int), color="orange", linestyle="--", linewidth=2, label="Ausgleichsgerade")

plt.xscale("log"); plt.yscale("log")
plt.xlabel("Frequenz [Hz]"); plt.ylabel(r"$V = \frac{U_{aus}}{U_{ein}}$")
#plt.title("Integrator: Amplitude vs Frequenz")
plt.legend(); plt.grid(True, linestyle="--", linewidth=0.5); plt.tight_layout()
plt.show()

# === Ergebnisse Integrator ausgeben ===
print("=== Integrator ===")
print(f"Nonlinearer Fit: a_int = {a_int:.6g} ± {a_err_int:.6g}, b_int = {b_int:.6g} ± {b_err_int:.6g}")
print(f"Ausgleichsgerade = {a_lin_int:.6g}, b_lin_int = {b_lin_int:.6g}")
print(f"Ckonst_int = {Ckonst_int:.6g} F")
print()

# -----------------------------
# === Differentiator ===
# -----------------------------
# === Eingangsspannung  ===
U_e_diff = 3.14  # Volt

# Parameter: C = 22 nF, R = 100 kΩ (wie du gesagt hast: 22 nF)
C_diff = 22e-9   # 22 nF
R_diff = 100e3   # 100 kΩ

# === Daten einlesen ===
f_diff, U_a_diff = np.genfromtxt("data/differentiator.txt", unpack=True)

# === Verstärkung berechnen ===
V_diff = U_a_diff / U_e_diff

# === Nonlinearer Fit mit curve_fit ===
p_diff, cov_diff = curve_fit(power_law7, f_diff, V_diff, maxfev=10000)
a_diff, b_diff = p_diff
a_err_diff, b_err_diff = np.sqrt(np.diag(cov_diff))

# === Lineare Ausgleichsgerade im ln-ln Raum ===
ln_f_diff = np.log(f_diff)
ln_V_diff = np.log(V_diff)
coef_lin_diff = np.polyfit(ln_f_diff, ln_V_diff, deg=1)
b_lin_diff = coef_lin_diff[0]
a_lin_diff = np.exp(coef_lin_diff[1])

# === Zeitkonstante berechnen (wie in deinem Ansatz) ===
Ckonst_diff = a_diff * R_diff * C_diff

# === Plot Differentiator: Daten, curve_fit (rot), ln-ln Ausgleich (orange gestrichelt) ===
plt.figure(figsize=(8,5))
plt.scatter(f_diff, V_diff, label="Differentiator Messdaten", color="blue", s=20)

f_plot_diff = np.logspace(np.log10(f_diff.min()*0.8), np.log10(f_diff.max()*1.2), 400)
plt.plot(f_plot_diff, power_law7(f_plot_diff, a_diff, b_diff), color="red", linewidth=2, label="Fit")
plt.plot(f_plot_diff, power_law7(f_plot_diff, a_lin_diff, b_lin_diff), color="orange", linestyle="--", linewidth=2, label="Ausgleichsgerade")

plt.xscale("log"); plt.yscale("log")
plt.xlabel("Frequenz [Hz]"); plt.ylabel(r"$V =\frac{U_{aus}}{U_{ein}}$")
#plt.title("Differentiator: Amplitude vs Frequenz")
plt.legend(); plt.grid(True, linestyle="--", linewidth=0.5); plt.tight_layout()
plt.show()

# === Ergebnisse Differentiator ausgeben ===
print("=== Differentiator ===")
print(f"Fit : a_diff = {a_diff:.6g} ± {a_err_diff:.6g}, b_diff = {b_diff:.6g} ± {b_err_diff:.6g}")
print(f"Ausgleichsgerade : a_lin_diff = {a_lin_diff:.6g}, b_lin_diff = {b_lin_diff:.6g}")
print(f"Ckonst_diff  = {Ckonst_diff:.6g} F")


