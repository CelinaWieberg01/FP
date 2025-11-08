import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
from scipy.signal import find_peaks


def pow(x, a, b):
    return a*x**b



def plot1(f_inv, U_inv, R1, R2, U_in_inv, c1, c2, i, error):
    U_inv /= U_in_inv
    const = c1
    powerlaw = c2

    f_const = f_inv[:const]
    U_const = U_inv[:const]
    err_const = error[:const]

    f_oth = f_inv[const:powerlaw]
    U_oth = U_inv[const:powerlaw]
    err_oth = error[const:powerlaw]

    f_pow = f_inv[powerlaw:]
    U_pow = U_inv[powerlaw:]
    err_pow = error[powerlaw:]

    plt.errorbar(f_const, U_const, yerr=err_const, fmt="o", capsize=2, color="green", label=r"Konst.")
    plt.errorbar(f_oth,   U_oth,   yerr=err_oth,   fmt="x", capsize=2, color="red",   label=r"ungefittet")
    plt.errorbar(f_pow,   U_pow,   yerr=err_pow,   fmt="o", capsize=2, color="blue",  label=r"Potenzgesetz")



    params, cov = curve_fit(pow, f_pow, U_pow)

    a = params[0]
    a_err = np.sqrt(np.diag(cov))[0]

    b = params[1]
    b_err = np.sqrt(np.diag(cov))[1]

    a_u = ufloat(a, a_err)
    b_u = ufloat(b, b_err)

    print(i, f"a = {a:.3} +- {a_err:.3}, b = {b:.3} +- {b_err:.3}")

    fit_x = np.linspace(f_pow[0]-f_pow[0]/2, f_pow[-1]+f_pow[-1]/2)

    plt.plot(fit_x, pow(fit_x, a, b), color="lightblue", label=r"Fit Potenzgesetz")

    # linke flanke
    plateau_height = np.mean(U_const)
    plateau_height_err = U_const.std(ddof=1)/np.sqrt(len(U_const))
    plateau_u = ufloat(plateau_height, plateau_height_err)

    plateau_x = np.linspace(f_const[0]-f_const[0]/2, f_const[-1]+f_const[-1]*10)

    print(i, f"Plateau = {plateau_height:.3} +- {plateau_height_err:.3}")

    plt.hlines(plateau_height, plateau_x[0], plateau_x[-1], color="lightgreen", label=r"Fit Potenzgesetz")

    plt.hlines(plateau_height/np.sqrt(2), f_inv[0], f_inv[-1], linestyles="dashed", color="orange", label=r"$\frac{U_\text{const}}{\sqrt{2}}$")

    f_cutoff = (plateau_u/np.sqrt(2) / a_u)**(1/b_u)
    print(i, "f_cutoff = ", f_cutoff)

    theoretical_amplification = R2/R1
    print(i, "A_theo = ", theoretical_amplification)

    dev_amplification = np.abs(theoretical_amplification - plateau_u)/theoretical_amplification
    print(i, "Deviation ampl = ", dev_amplification)
    plt.xscale("log")
    plt.xlabel(r"$f$ in $\si{\hertz}$")

    plt.yscale("log")
    plt.ylabel(r"$A$")

    #plt.title(r"R_1 = \SI{}, R2 = {R2}")
    plt.legend(loc="lower left")
    plt.grid(True, which="both")

make_plots = False

if make_plots == True:
    print("\n\n\n")
    # data 1
    f_inv, U_inv = np.genfromtxt("data/Invamp1.txt", unpack=True) # Hz, V
    R1 = 1000
    R2 = 100000

    U_in_inv = 107e-3 # V
    U_in_inv_u = ufloat(U_in_inv, 10e-3)
    U_inv_u = unp.uarray(U_inv, 0.05)

    errs = unp.std_devs(U_inv_u/U_in_inv_u)

    c1 = 7
    c2 = 13

    plot1(f_inv, U_inv, R1, R2, U_in_inv, c1, c2, 1, errs)
    plt.title(r"Inverting Amplifier mit $R_1 = \SI{1}{\kilo\ohm}$ und $R_2 = \SI{100}{\kilo\ohm}$")
    plt.savefig("plots/inv1.pdf")
    plt.figure()

    print("\n\n\n")
    # data 2
    f_inv, U_inv = np.genfromtxt("data/Invamp2.txt", unpack=True) # Hz, V
    R1 = 1000
    R2 = 10000
    U_in_inv = 107e-3 # V
    U_in_inv_u = ufloat(U_in_inv, 10e-3)
    U_inv_u = unp.uarray(U_inv, 0.05)

    errs = unp.std_devs(U_inv_u/U_in_inv_u)

    c1 = 5
    c2 = -6

    plot1(f_inv, U_inv, R1, R2, U_in_inv, c1, c2, 2, errs)
    plt.title(r"Inverting Amplifier mit $R_1 = \SI{1}{\kilo\ohm}$ und $R_2 = \SI{10}{\kilo\ohm}$")
    plt.savefig("plots/inv2.pdf")
    plt.figure()

    print("\n\n\n")
    # data 3
    f_inv, U_inv = np.genfromtxt("data/Invamp3.txt", unpack=True) # Hz, V
    R1 = 1000
    R2 = 150000
    U_in_inv = 40e-3 # V
    U_in_inv_u = ufloat(U_in_inv, 10e-3)
    U_inv_u = unp.uarray(U_inv, 0.05)

    errs = unp.std_devs(U_inv_u/U_in_inv_u)

    c1 = 6
    c2 = -8

    plot1(f_inv, U_inv, R1, R2, U_in_inv, c1, c2, 3, errs)
    plt.title(r"Inverting Amplifier mit $R_1 = \SI{1}{\kilo\ohm}$ und $R_2 = \SI{150}{\kilo\ohm}$")
    plt.savefig("plots/inv3.pdf")
    plt.figure()




    print("\n\n\n")
    # integrator

    def lin(x, m, b):
        return m*x+b

    f_int, U1, U2 = np.genfromtxt("data/IntegreatorAmp.txt", unpack=True)
    U1 *= 1e-3
    U1 = unp.uarray(U1, 0.005)
    U2 = unp.uarray(U2, 0.05)
    U_int = U2/U1

    x = 1/(2*np.pi*f_int)

    R = 10000
    C = 100e-9
    tau = R*C
    print("INT: theoretische Zeitkonstante tau = ", tau)

    params, cov = curve_fit(lin, x, unp.nominal_values(U_int), sigma=unp.std_devs(U_int))

    m = params[0]
    m_err = np.sqrt(np.diag(cov))[0]

    b = params[1]
    b_err = np.sqrt(np.diag(cov))[1]

    m_u = ufloat(m, m_err)
    b_u = ufloat(b, b_err)

    print(f"INT: m = {m:.3} +- {m_err:.3}, b = {b:.3} +- {b_err:.3}")

    tau_fit = 1/m_u
    print("INT: Zeitkonstante im Fit tau = ", tau_fit)

    dev_tau = np.abs(tau - tau_fit)/tau
    print("INT: deviation tau = ", dev_tau)

    int_x = np.linspace(x[0], x[-1])
    plt.plot(int_x, lin(int_x, *params), color="lightblue", label=r"Fit")



    plt.errorbar(x, unp.nominal_values(U_int), yerr=unp.std_devs(U_int), fmt="o", capsize=2, color="blue", label=r"Messung")

    #plt.xscale("log")
    plt.xlabel(r"$\frac{1}{\omega} = \frac{1}{2\pi f}$ in $\si{\per\hertz}$")

    #plt.yscale("log")
    plt.ylabel(r"A")
    plt.legend()
    plt.title("Integrator")
    plt.grid("on")
    plt.savefig("plots/int.pdf")
    plt.figure()




    print("\n\n\n")
    # differentiatior

    f_int, U2 = np.genfromtxt("data/differentiator.txt", unpack=True)
    U1 = 3.14
    U1_u = ufloat(U1, 0.05)
    U2 = unp.uarray(U2, 0.05)
    U_int = U2/U1

    x = 2*np.pi*f_int

    R = 100000
    C = 22e-9
    tau = R*C
    print("DIFF: theoretische Zeitkonstante tau = ", tau)

    params, cov = curve_fit(lin, x, unp.nominal_values(U_int), sigma=unp.std_devs(U_int))

    m = params[0]
    m_err = np.sqrt(np.diag(cov))[0]

    b = params[1]
    b_err = np.sqrt(np.diag(cov))[1]

    m_u = ufloat(m, m_err)
    b_u = ufloat(b, b_err)

    print(f"DIFF: m = {m:.3} +- {m_err:.3}, b = {b:.3} +- {b_err:.3}")

    tau_fit = m_u
    print("DIFF: Zeitkonstante im Fit tau = ", tau_fit)

    dev_tau = np.abs(tau - tau_fit)/tau
    print("DIFF: deviation tau = ", dev_tau)

    int_x = np.linspace(x[0], x[-1])
    plt.plot(int_x, lin(int_x, *params), color="lightblue", label=r"Fit")



    plt.errorbar(x, unp.nominal_values(U_int), yerr=unp.std_devs(U_int), fmt="o", capsize=2, color="blue", label=r"Messung")

    #plt.xscale("log")
    plt.xlabel(r"$\omega = 2\pi f$ in $\si{\hertz}$")

    #plt.yscale("log")
    plt.ylabel(r"$A$")
    plt.legend()
    plt.grid("on")
    plt.title(r"Differentiator")
    plt.savefig("plots/diff.pdf")
    plt.figure()




    print("\n\n\n")
    # Schmitt Trigger
    R1 = 10000
    R2 = 100000
    U_s = ufloat(28.1, 0.1)
    U_pm = R1/R2 * U_s
    print("SCHMITT: Soll-Schwelle = ", U_pm)

    U_meas = ufloat(2.980, 0.001)

    print("SCHMITT: Messung = ", U_meas)

    dev = np.abs(U_pm - U_meas)/U_pm
    print("SCHMITT: Dev = ", dev)



    print("\n\n\n")
    # generator 1
    R1 = 10000
    R2 = 100000
    R3 = 1000
    C = 1e-6

    f = R2 / (4 * C * R1 * R3)
    print("GEN 1: Theoretische Frequenz = ", f)

    f_meas = ufloat(1645.6, 100)
    print("GEN 1: gemessene Frequenz = ", f_meas)

    dev = np.abs(f - f_meas)/f
    print("GEN 1: Dev = ", dev)



    print("\n\n\n")
    # generator 2
    C = 100e-9 
    R = 10000

    T = 2*np.pi*R*C
    tau = 20*R*C

    print(f"GEN 2: theoretisches T = {T:.3}")
    print(f"GEN 2: theoretisches tau = {tau:.3}")



    def decay(x, a, tau):
        return a*np.exp(-x/tau)


    t, U = np.genfromtxt("data/csv/scope_10.csv", delimiter=",", usecols=(0,2), unpack=True)

    plt.plot(t, U, color="blue", linewidth=1, label=r"Messung")
    plt.xlabel(r"$t$ in $\si{\second}$")
    plt.ylabel(r"rel. Spannung $U$ in $\si{\volt}$")
    plt.legend()
    plt.grid("on")
    plt.title(r"Generator mit Dämpfung bei maximalem Potentiometer")
    plt.savefig("plots/gen2_max.pdf")
    plt.figure()

    t = t[U > 0.03]
    U = U[U > 0.03]

    U = U[t < 0.3]
    t = t[t < 0.3]

    U -= np.mean(U)
    t -= t[0]

    pos_peaks, _ = find_peaks(U, distance=15)
    t_ppeaks = t[pos_peaks]
    U_ppeaks = U[pos_peaks]
    t_ppeaks_u = unp.uarray(t_ppeaks, 0.0015)
    U_ppeaks_u = unp.uarray(U_ppeaks, 0.001)


    neg_peaks, _ = find_peaks(-U, distance=15)
    t_npeaks = t[neg_peaks]
    U_npeaks = U[neg_peaks]
    t_npeaks_u = unp.uarray(t_npeaks, 0.0015)
    U_npeaks_u = unp.uarray(U_npeaks, 0.0015)

    plt.errorbar(t_ppeaks, U_ppeaks, yerr=0.0015, xerr=0.0015, fmt="x", capsize=2, color="red", label=r"pos. Peaks")
    plt.errorbar(t_npeaks, U_npeaks, yerr=0.0015, xerr=0.0015, fmt="x", capsize=2, color="green", label=r"pos. Peaks")

    Tp = t_ppeaks[1:] - t_ppeaks[:-1]
    Tn = t_npeaks[1:] - t_npeaks[:-1]
    Ts = np.concatenate((Tp, Tn))
    T_mean = np.mean(Ts)
    T_err = np.std(Ts)/np.sqrt(len(Ts))
    T_u = ufloat(T_mean, T_err)
    print("GEN 2_max: T = ", T_u)

    t_all_peaks = np.concatenate((t_ppeaks, t_npeaks))
    U_all_peaks = np.concatenate((U_ppeaks, -U_npeaks))
    params, cov = curve_fit(decay, t_all_peaks, U_all_peaks, sigma=0.0015, p0=(-0.01, 0.10))
    errs = np.sqrt(np.diag(cov))
    print(f"GEN 2_max: a_n = {params[0]:.3} +- {errs[0]:.3}, tau_n = {params[1]:.3} +- {errs[1]:.3}")

    plt.plot(t, U, label=r"Messung")
    plt.xlabel(r"$t$ in $\si{\second}$")
    plt.ylabel(r"rel. Spannung $U$ in $\si{\volt}$")
    plt.title(r"Generator mit Dämpfung bei maximalem Potentiometer, Ausschnitt")
    plt.legend()
    plt.grid("on")
    plt.savefig("plots/gen2_max_sec.pdf")
    plt.figure()

    dev_T = np.abs(T-T_u)/T
    print(f"GEN 2_max: dev_T = {dev_T}")
    dev_tau = np.abs(tau - ufloat(params[1], errs[1]))/tau
    print(f"GEN 2_max: dev_tau = {dev_tau}")

    print("\n\n\n")
    t, U = np.genfromtxt("data/csv/scope_13.csv", delimiter=",", usecols=(0,2), unpack=True)

    plt.plot(t, U, color="blue", linewidth=1, label=r"Messung")
    plt.xlabel(r"$t$ in $\si{\second}$")
    plt.ylabel(r"rel. Spannung $U$ in $\si{\volt}$")
    plt.legend()
    plt.grid("on")
    plt.title(r"Generator mit Dämpfung bei minimalem Potentiometer")
    plt.savefig("plots/gen2_min.pdf")
    plt.figure()

    t = t[U > 0.03]
    U = U[U > 0.03]

    U = U[t < 0.28]
    t = t[t < 0.28]

    U -= np.mean(U)
    t -= t[0]


    pos_peaks, _ = find_peaks(U, distance=15)
    t_ppeaks = t[pos_peaks]
    U_ppeaks = U[pos_peaks]
    t_ppeaks_u = unp.uarray(t_ppeaks, 0.0015)
    U_ppeaks_u = unp.uarray(U_ppeaks, 0.001)


    neg_peaks, _ = find_peaks(-U, distance=15)
    t_npeaks = t[neg_peaks]
    U_npeaks = U[neg_peaks]
    t_npeaks_u = unp.uarray(t_npeaks, 0.0015)
    U_npeaks_u = unp.uarray(U_npeaks, 0.0015)

    plt.errorbar(t_ppeaks, U_ppeaks, yerr=0.0015, xerr=0.0015, fmt="x", capsize=2, color="red", label=r"pos. Peak")
    plt.errorbar(t_npeaks, U_npeaks, yerr=0.0015, xerr=0.0015, fmt="x", capsize=2, color="green", label=r"pos. Peak")

    Tp = t_ppeaks[1:] - t_ppeaks[:-1]
    Tn = t_npeaks[1:] - t_npeaks[:-1]
    Ts = np.concatenate((Tp, Tn))
    T_mean = np.mean(Ts)
    T_err = np.std(Ts)/np.sqrt(len(Ts))
    T_u = ufloat(T_mean, T_err)
    print("GEN 2_min: T = ", T_u)

    t_all_peaks = np.concatenate((t_ppeaks, t_npeaks))
    U_all_peaks = np.concatenate((U_ppeaks, -U_npeaks))
    params, cov = curve_fit(decay, t_all_peaks, U_all_peaks, sigma=0.0015, p0=(0.01, 0.01))
    errs = np.sqrt(np.diag(cov))
    print(f"GEN 2_min: a_n = {params[0]:.3} +- {errs[0]:.3}, tau_n = {params[1]:.3} +- {errs[1]:.3}")

    plt.plot(t, U, label=r"Messung")
    plt.xlabel(r"$t$ in $\si{\second}$")
    plt.ylabel(r"rel. Spannung $U$ in $\si{\volt}$")
    plt.title(r"Generator mit Dämpfung bei minimalem Potentiometer, Ausschnitt")
    plt.legend()
    plt.grid("on")
    plt.savefig("plots/gen2_min_sec.pdf")

    dev_T = np.abs(T-T_u)/T
    print(f"GEN 2_min: dev_T = {dev_T}")
    dev_tau = np.abs(tau - ufloat(params[1], errs[1]))/tau
    print(f"GEN 2_min: dev_tau = {dev_tau}")