import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

DP = False
TP = False
# Funktion für Linfits
def linfit(x, a, b):
    return a*x + b

#################################
    #### Drehschieberpumpe

if DP == True:
    print("DP läuft")
    # t und D0 = usecols(2, 3)

    # Fehler für Evakuierung
    # 1000 - 10 mbar: 0.3%
    # 10 - 2*10^-3:    10%
    pE = ufloat(3.83E-3, 3.83E-3 * 0.1) # Enddruck
    VD = ufloat(34, 34*0.1) # Volumen Drehschieberpumpe

    t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Evak_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0] # Verschieb Zeit 
    D0_err = np.concatenate((D0[D0>=10]*0.003, D0[D0<10]*0.1)) # D0_err für Plots (nicht nützlich weil kaum erkennbar)

    #print("DP Evakuierungskurve")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 20 == 0:
    #        print(f"{time} & {D0[index]} & {D0_err[index]}")
    #print(" ")
    """
DP Evakuierungskurve
t     p(t)     err
0 & 1020.0 & 3.06
20 & 441.8 & 1.32
40 & 208.7 & 0.63
60 & 91.3 & 0.27
80 & 39.5 & 0.11
100 & 17.5 & 0.05
120 & 7.7 & 0.77
140 & 3.5 & 0.35
160 & 2.0 & 0.20
180 & 1.3 & 0.13
200 & 0.92 & 0.09
220 & 0.69 & 0.07
240 & 0.57 & 0.06
260 & 0.46 & 0.05
280 & 0.39 & 0.04
300 & 0.33 & 0.03
320 & 0.29 & 0.03
340 & 0.25 & 0.025
360 & 0.22 & 0.02
380 & 0.20 & 0.02
400 & 0.18 & 0.02
420 & 0.16 & 0.02
440 & 0.14 & 0.01
460 & 0.13 & 0.01
480 & 0.12 & 0.01
500 & 0.11 & 0.011
520 & 0.10 & 0.01
540 & 0.094 & 0.0094
560 & 0.087 & 0.0087
580 & 0.079 & 0.0079
600 & 0.073 & 0.0073
620 & 0.067 & 0.0067
    """

    ### Allgemeiner Plot 

    plt.errorbar(t, D0, yerr=D0_err, errorevery=(10), fmt=".", ecolor="lightblue", label="Messung")

    # Zeiten für LinFit: 0 - 150, 150 - 250, 250 - Ende 
    # 0 - 150
    params, cov = curve_fit(linfit, t[0:150], np.log((D0[0:150]-pE)/(D0[0]-pE)), sigma=D0_err[0:150], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x0150 = np.linspace(0, 170)
    plt.errorbar(x0150, pE + (D0[0] - pE)*np.exp(linfit(x0150, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
    #print("DP Fit 1 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    # 150 - 250
    params, cov = curve_fit(linfit, t[150:250], np.log((D0[150:250]-pE)/(D0[0]-pE)), sigma=D0_err[150:250], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x150250 = np.linspace(130, 300)
    plt.errorbar(x150250, pE + (D0[0] - pE)*np.exp(linfit(x150250, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
    #print("DP Fit 2 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    # 250 - Ende
    params, cov = curve_fit(linfit, t[250:], np.log((D0[250:]-pE)/(D0[0]-pE)), sigma=D0_err[250:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VD
    x250end = np.linspace(230, t[-1])
    plt.errorbar(x250end, pE + (D0[0] - pE)*np.exp(linfit(x250end, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
    #print("DP Fit 3 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    plt.yscale("log")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Evakuierungskurve Drehschieberpumpe")
    plt.legend()

    plt.savefig("DP_Evakuierungskurve.pdf")
    plt.figure()
    #print(" ")

    ################################

    # Fehler für Leckratenmessung
    # 0.5 mbar:        10%
    # Sonst:          0.3%

    # 0.5 mbar
    pG = unp.uarray(0.5, 0.05)
    t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_05mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D0_err = D0*0.1

    plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VD/pG * linfit_params[0]
    #print("DP Leck 0.5 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    #print("DP Leckrate 5 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D0[index]} & {D0_err[index]}")
    #print(" ")

    """
DP Leckrate 5 mbar
t     p(t)     err
0 & 0.51 & 0.05
10 & 1.8 & 0.18
20 & 2.0 & 0.20
30 & 2.1 & 0.21
40 & 2.2 & 0.22
50 & 2.4 & 0.24
60 & 2.5 & 0.25
70 & 2.7 & 0.27
80 & 2.8 & 0.28
90 & 3.0 & 0.30
100 & 3.1 & 0.31
110 & 3.2 & 0.32
120 & 3.4 & 0.34
130 & 3.6 & 0.36
140 & 3.8 & 0.38
150 & 3.9 & 0.39
160 & 4.1 & 0.41
170 & 4.2 & 0.42
180 & 4.3 & 0.43
190 & 4.5 & 0.45
200 & 4.6 & 0.46
210 & 4.8 & 0.48
    """

    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $0,5$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("DP_Leck_05mbar.pdf")
    plt.figure()
###############################################################

    # 10 mbar
    pG = unp.uarray(10, 10*0.03)
    t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_10mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D0_err = D0*0.003

    plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VD/pG * linfit_params[0]
    #print("DP Leck 10 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    #print("DP Evakuierungskurve 10 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D0[index]} & {D0_err[index]}")
    #print(" ")
    """
    DP Leckrate 10 mbar
t     p(t)     err
0 & 10.0 & 0.03
10 & 20.0 & 0.06
20 & 24.7 & 0.0741
30 & 28.2 & 0.0846
40 & 33.7 & 0.1011
50 & 38.0 & 0.114
60 & 43.1 & 0.1293
70 & 47.7 & 0.1431
80 & 52.4 & 0.1572
90 & 57.1 & 0.1713
100 & 61.4 & 0.1842
110 & 66.4 & 0.1992
120 & 71.1 & 0.2133
130 & 76.1 & 0.2283
140 & 81.1 & 0.2433
150 & 85.8 & 0.2574
160 & 90.5 & 0.2715
170 & 95.6 & 0.2868
180 & 100.6 & 0.3018
190 & 105.7 & 0.3171
200 & 110.7 & 0.3321
    """

    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $10$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("DP_Leck_10mbar.pdf")
    plt.figure()
################################################

    # 50 mbar
    pG = unp.uarray(50, 50*0.03)
    t, D0 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck_50mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D0_err = D0*0.003

    plt.errorbar(t, D0, yerr=D0_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    params, cov = curve_fit(linfit, t[10:], D0[10:], sigma=D0_err[10:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VD/pG * linfit_params[0]
    #print("DP Leck 50 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    #print("DP Leckrate 50 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D0[index]} & {D0_err[index]}")
    #print(" ")

    """
    DP Leckrate 50 mbar
t     p(t)     err
0 & 51.2 & 0.1536
10 & 85.4 & 0.2562
20 & 109.4 & 0.3282
30 & 131.3 & 0.3939
40 & 155.1 & 0.4653
50 & 178.9 & 0.5367
60 & 203.7 & 0.6111
70 & 227.5 & 0.6825
80 & 253.2 & 0.7596
90 & 275.1 & 0.8253
100 & 299.0 & 0.897
110 & 324.6 & 0.9738
120 & 348.5 & 1.0455
130 & 374.2 & 1.1226
140 & 399.7 & 1.1991
150 & 425.2 & 1.2756
160 & 448.8 & 1.3464
170 & 474.0 & 1.422
180 & 499.2 & 1.4976
190 & 524.1 & 1.5723
200 & 548.9 & 1.6467
    """

    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $50$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("DP_Leck_50mbar.pdf")
    plt.figure()

    # 100 mbar 1
    pG = unp.uarray(100, 100*0.03)
    t1, D01 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck1_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t1 = t1 - t1[0]
    D01_err = D01*0.003
    D01 = unp.uarray(D01, D01_err)

    t2, D02 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck2_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t2 = t2 - t2[0]
    D02_err = D02*0.003
    D02 = unp.uarray(D02, D02_err)
    
    t3, D03 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Leck3_100mbar_D0.csv", delimiter=",", skip_header=1, usecols=(2,3), dtype=None, encoding=None, unpack=True) # Unpack Data
    t3 = t3 - t3[0]
    D03_err = D03*0.003
    D03 = unp.uarray(D03, D03_err)

    plt.errorbar(t1, unp.nominal_values(D01), yerr=D01_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung 1")
    plt.errorbar(t2, unp.nominal_values(D02), yerr=D02_err, fmt=".", ecolor="yellow", errorevery=(5), label="Messung 2")
    plt.errorbar(t3, unp.nominal_values(D03), yerr=D03_err, fmt=".", ecolor="lightgreen", errorevery=(5), label="Messung 3")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $100$ $\mathrm{mbar}$, alle Messungen")
    plt.legend()
    plt.savefig("DP_Leck_100mbar_alle.pdf")
    plt.figure()

    # Mittelwertbildung 100 mbar
    N = min(len(t1), len(t2), len(t3))
    D01 = D01[:N]
    D02 = D02[:N]
    D03 = D03[:N]
    mean = (unp.nominal_values(D01) + unp.nominal_values(D02) + unp.nominal_values(D03))/3
    std_all = np.array((unp.std_devs(D01), unp.std_devs(D02), unp.std_devs(D03)))
    std_mean = np.sqrt(np.sum(std_all**2, axis=0)/N)
    t = t1[:N]
    D0 = unp.uarray(mean, std_mean)

        # plot für 100 mbar mittelwerte
    plt.errorbar(t, unp.nominal_values(D0), yerr=unp.std_devs(D0), fmt=".", ecolor="lightblue", errorevery=(5), label="Mittelwert")

        # fit für 100 mbar mittelwerte
    params, cov = curve_fit(linfit, t[10:150], unp.nominal_values(D0[10:150]), sigma=unp.std_devs(D0[10:150]), absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VD/pG * linfit_params[0]
    #print("DP Leck Mittel, 100 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    
    #print("DP Leckrate 100 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate((t3)):
    #    if time % 10 == 0:
    #        print(f"{time} & {unp.nominal_values(D01[index])} & {D01_err[index]} & {unp.nominal_values(D02[index])} & {D02_err[index]} & {unp.nominal_values(D03[index])} & {D03_err[index]} & {mean[index]} & {std_mean[index]}")
    #print(" ")

    """
    DP Leckrate 100 mbar
t     p1     err    p2   err    p3      err     p_mean      err
0 & 101.9 & 0.3057 & 99.1 & 0.2973 & 99.5 & 0.2985 & 100.16666666666667 & 0.03653340667054833
10 & 165.9 & 0.49770000000000003 & 156.5 & 0.46950000000000003 & 155.7 & 0.46709999999999996 & 159.36666666666665 & 0.058145305923911185
20 & 210.9 & 0.6327 & 201.1 & 0.6033 & 202.8 & 0.6084 & 204.9333333333333 & 0.0747551253659099
30 & 253.6 & 0.7608 & 240.8 & 0.7224 & 240.8 & 0.7224 & 245.0666666666667 & 0.0894025093045533
40 & 300.0 & 0.9 & 285.3 & 0.8559000000000001 & 289.3 & 0.8679 & 291.5333333333333 & 0.10634582540598755
50 & 349.8 & 1.0494 & 333.3 & 0.9999 & 334.4 & 1.0031999999999999 & 339.1666666666667 & 0.12372404161239804
60 & 396.1 & 1.1883000000000001 & 377.7 & 1.1331 & 379.3 & 1.1379000000000001 & 384.3666666666666 & 0.14021076517140452
70 & 445.4 & 1.3362 & 425.4 & 1.2762 & 430.8 & 1.2924 & 433.8666666666666 & 0.15826047370378307
80 & 494.3 & 1.4829 & 469.4 & 1.4082 & 481.8 & 1.4454 & 481.8333333333333 & 0.1757629506681832
90 & 535.5 & 1.6065 & 509.3 & 1.5279 & 521.9 & 1.5656999999999999 & 522.2333333333332 & 0.19049762614400859
100 & 579.6 & 1.7388000000000001 & 555.4 & 1.6662 & 568.1 & 1.7043000000000001 & 567.6999999999999 & 0.20707065939020808
110 & 625.9 & 1.8777 & 600.2 & 1.8006000000000002 & 613.3 & 1.8398999999999999 & 613.1333333333333 & 0.22364148892033084
120 & 670.7 & 2.0121 & 640.9 & 1.9227 & 657.2 & 1.9716000000000002 & 656.2666666666667 & 0.2393806485011838
130 & 716.4 & 2.1492 & 683.1 & 2.0493 & 699.4 & 2.0982 & 699.6333333333333 & 0.25520333214696916
140 & 756.7 & 2.2701000000000002 & 723.5 & 2.1705 & 739.4 & 2.2182 & 739.8666666666667 & 0.26987349253019155
150 & 794.6 & 2.3838 & 764.3 & 2.2929 & 777.2 & 2.3316000000000003 & 778.7000000000002 & 0.28402674031122965
160 & 827.3 & 2.4819 & 794.9 & 2.3847 & 809.9 & 2.4297 & 810.6999999999999 & 0.29570041914372286
170 & 859.6 & 2.5788 & 828.4 & 2.4852 & 844.7 & 2.5341 & 844.2333333333332 & 0.30792561389819667
180 & 890.7 & 2.6721000000000004 & 861.2 & 2.5836 & 874.5 & 2.6235 & 875.4666666666667 & 0.3193116026797513
190 & 918.1 & 2.7543 & 888.9 & 2.6667 & 903.1 & 2.7093000000000003 & 903.3666666666667 & 0.3294850727312383
200 & 940.2 & 2.8206 & 913.7 & 2.7411000000000003 & 928.2 & 2.7846 & 927.3666666666668 & 0.33823222507445644
    """

    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $100$ $\mathrm{mbar}$, Mittelwert")
    plt.legend()
    plt.savefig("DP_Leck1_100mbar_mittelwert.pdf")
    plt.figure()



'''
DP Fit 1 S =  1.39+/-0.14  l/s
m =  -0.0408+/-0.0007    b =  0.01+/-0.07
 
DP Fit 2 S =  0.49+/-0.05  l/s
m =  -0.01441+/-0.00034    b =  -4.16+/-0.08
 
DP Fit 3 S =  0.41+/-0.04  l/s
m =  -0.011943+/-0.000008    b =  -4.116+/-0.004
 
DP Fit 4 S =  0.203+/-0.020  l/s
m =  -0.005971+/-0.000004    b =  -6.174+/-0.006
 
 
DP Leck 0.5 mbar S =  1.02+/-0.15  l/s
m =  0.0150+/-0.0004    b =  1.63+/-0.04
 
DP Leck 10 mbar S =  1.60+/-0.17  l/s
m =  0.47184+/-0.00022    b =  14.866+/-0.016
 
DP Leck 50 mbar S =  1.64+/-0.17  l/s
m =  2.4190+/-0.0010    b =  59.78+/-0.07
 
DP Leck Mittel, 100 mbar S =  1.54+/-0.16  l/s
m =  4.52830+/-0.00032    b =  112.640+/-0.018
'''
####################################################
####################################################
####################################################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
# 100 - 10^-8 mbar: 30%

if TP == True:
    print("TP läuft")
    pE = 5.0E-6 # Enddruck
    VT = ufloat(33, 33*0.1) # Volumen Drehschieberpumpe

    t1, D21 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak1_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t1 = t1 - t1[0] # Verschieb Zeit 
    D21_err = D21*0.3

    t2, D22 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak2_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t2 = t2 - t2[0] # Verschieb Zeit 
    D22_err = D22*0.3

    t3, D23 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Evak3_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t3 = t3 - t3[0] # Verschieb Zeit 
    D23_err = D23*0.3

    ### Allgemeiner Plot 

    plt.errorbar(t1, D21, yerr=D21_err, errorevery=(5), fmt=".", ecolor="lightblue", label="Messung 1")
    plt.errorbar(t2, D22, yerr=D22_err, errorevery=(5), fmt=".", ecolor="yellow", label="Messung 2")
    plt.errorbar(t3, D23, yerr=D23_err, errorevery=(5), fmt=".", ecolor="lightgreen", label="Messung 3")
    plt.yscale("log")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Evakuierungskurve Turbomolekularpumpe, alle Messungen")
    plt.legend()
    plt.savefig("TP_Evakuierungskurve_alle.pdf")
    plt.figure()

    # Mittelwertbildung TP
    N = min(len(t2), len(t3))
    D22 = D22[:N]
    D23 = D23[:N]

    D22_err = D22_err[:N]
    D23_err = D23_err[:N]

    mean = (D22 + D23)/3
    std_all = np.array((D22_err, D23_err))
    std_mean = np.sqrt(np.sum(std_all**2, axis=0)/N)
    t = t1[:N]
    D2 = mean
    D2_err = std_mean

    #print("TP Evak")
    #print("t     p1     err     p2      err     p3      err     p_")
    #for index, time in enumerate((t3)):
    #    if time % 10 == 0:
    #        print(f"{time} & {unp.nominal_values(D01[index])} & {D01_err[index]} & {unp.nominal_values(D02[index])} & {D02_err[index]} & {unp.nominal_values(D03[index])} & {D03_err[index]} & {mean[index-2]} & {std_mean[index-2]}")
    #print(" ")
    """
    TP Evak
t     p1     err     p2      err     p3      err     p_
0 & 101.9 & 0.3057 & 99.1 & 0.2973 & 99.5 & 0.2985 & 8.6e-06 & 4.857091576517491e-07
10 & 165.9 & 0.49770000000000003 & 156.5 & 0.46950000000000003 & 155.7 & 0.46709999999999996 & 7.999999999999999e-05 & 4.527709959826817e-06
20 & 210.9 & 0.6327 & 201.1 & 0.6033 & 202.8 & 0.6084 & 1.826666666666667e-05 & 1.0319769270447647e-06
30 & 253.6 & 0.7608 & 240.8 & 0.7224 & 240.8 & 0.7224 & 1.37e-05 & 7.74168082829847e-07
40 & 300.0 & 0.9 & 285.3 & 0.8559000000000001 & 289.3 & 0.8679 & 1.2066666666666666e-05 & 6.81791295103377e-07
50 & 349.8 & 1.0494 & 333.3 & 0.9999 & 334.4 & 1.0031999999999999 & 1.1066666666666665e-05 & 6.252294067169442e-07
60 & 396.1 & 1.1883000000000001 & 377.7 & 1.1331 & 379.3 & 1.1379000000000001 & 1.04e-05 & 5.874916652114422e-07
70 & 445.4 & 1.3362 & 425.4 & 1.2762 & 430.8 & 1.2924 & 9.933333333333334e-06 & 5.611475868015724e-07
80 & 494.3 & 1.4829 & 469.4 & 1.4082 & 481.8 & 1.4454 & 9.566666666666666e-06 & 5.404001142344133e-07
100 & 579.6 & 1.7388000000000001 & 555.4 & 1.6662 & 568.1 & 1.7043000000000001 & 9.1e-06 & 5.139724861657801e-07
110 & 625.9 & 1.8777 & 600.2 & 1.8006000000000002 & 613.3 & 1.8398999999999999 & 8.9e-06 & 5.026802180062421e-07
120 & 670.7 & 2.0121 & 640.9 & 1.9227 & 657.2 & 1.9716000000000002 & 8.766666666666667e-06 & 4.951520883033627e-07
130 & 716.4 & 2.1492 & 683.1 & 2.0493 & 699.4 & 2.0982 & 8.6e-06 & 4.857091576517491e-07
    """
    # plot für TP mittelwerte
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Mittelwert")

    # Zeiten für LinFit: 0 - 15, 15 - 50, 50 - Ende 
    # 0 - 15

    params, cov = curve_fit(linfit, t[:15], np.log((D2[:15]-pE)/(D2[0]-pE)), sigma=D2_err[:15], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VT
    x015 = np.linspace(0, 20)
    plt.errorbar(x015, pE + (D2[0] - pE)*np.exp(linfit(x015, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="red", label="Fit 1",)
    #print("TP Fit 1 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    # 15 - 50
    params, cov = curve_fit(linfit, t[15:50], np.log((D2[15:50]-pE)/(D2[0]-pE)), sigma=D2_err[15:50], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VT
    x1550 = np.linspace(10, 60)
    plt.errorbar(x1550, pE + (D2[0] - pE)*np.exp(linfit(x1550, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="green", label="Fit 2",)
    #print("TP Fit 2 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    # 50 - Ende
    params, cov = curve_fit(linfit, t[50:], np.log((D2[50:]-pE)/(D2[0]-pE)), sigma=D2_err[50:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    S = -linfit_params[0]*VT
    x50 = np.linspace(40, t[-1])
    plt.errorbar(x50, pE + (D2[0] - pE)*np.exp(linfit(x50, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1]))), color="purple", label="Fit 3",)
    #print("TP Fit 3 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    plt.yscale("log")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Evakuierungskurve Turbomolekularpumpe, Mittelwert")
    plt.legend()
    plt.savefig("TP_Evakuierungskurve.pdf")
    plt.figure()
######################################################################################
    # Fehler für Leckratenmessung
    # 100 - 10^-8 mbar: 30%
    # 
    # # 1e-4 mbar
    pG = unp.uarray(1E-4, (1E-4)*0.3)
    t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_1e4_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D2_err = D2*0.3
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    #print("TP Leckrate 1e-4 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D2[index]} & {D2_err[index]}")
    #print(" ")
    """
    TP Leckrate 1e-4 mbar
t     p(t)     err
0 & 0.000105 & 3.15e-05
10 & 0.00026 & 7.799999999999999e-05
20 & 0.000426 & 0.0001278
30 & 0.0006 & 0.00017999999999999998
40 & 0.0009 & 0.00027
50 & 0.00123 & 0.00036899999999999997
60 & 0.00159 & 0.000477
70 & 0.00203 & 0.0006090000000000001
80 & 0.00245 & 0.000735
90 & 0.00296 & 0.0008879999999999999
100 & 0.00346 & 0.0010379999999999999
110 & 0.00402 & 0.001206
120 & 0.0046 & 0.00138
    """

    params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VT/pG * linfit_params[0]
    #print("TP Leck 1E-4 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $10^{-4}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("TP_Leck_1e4.pdf")
    plt.figure()

    # # 2e-4 mbar
    pG = unp.uarray(2E-4, (2E-4)*0.3)
    t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_2e4_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D2_err = D2*0.3
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    #print("TP Leckrate 2e-4 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D2[index]} & {D2_err[index]}")
    #print(" ")
    """
    TP Leckrate 2e-4 mbar
t     p(t)     err
0 & 0.000206 & 6.18e-05
10 & 0.000678 & 0.00020339999999999998
20 & 0.00159 & 0.000477
30 & 0.00275 & 0.0008249999999999999
40 & 0.00413 & 0.0012389999999999999
50 & 0.00573 & 0.0017189999999999998
60 & 0.00715 & 0.0021449999999999998
70 & 0.00886 & 0.0026579999999999998
80 & 0.0105 & 0.00315
90 & 0.0116 & 0.0034799999999999996
100 & 0.0125 & 0.00375
110 & 0.0137 & 0.00411
120 & 0.0149 & 0.00447
    """

    params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VT/pG * linfit_params[0]
    #print("TP Leck 2E-4 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $2 \cdot 10^{-4}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("TP_Leck_2e4.pdf")
    plt.figure()

    # # 5e-5 mbar
    pG = unp.uarray(5E-5, (5E-5)*0.3)
    t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_5e5_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D2_err = D2*0.3
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    #print("TP Leckrate 5e-5 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D2[index]} & {D2_err[index]}")
    #print(" ")
    """
TP Leckrate 5e-5 mbar
    t     p(t)     err
0 & 5.01e-05 & 1.5029999999999998e-05
10 & 0.000124 & 3.72e-05
20 & 0.000175 & 5.2499999999999995e-05
30 & 0.000252 & 7.56e-05
40 & 0.000326 & 9.78e-05
50 & 0.000402 & 0.0001206
60 & 0.000469 & 0.0001407
70 & 0.000553 & 0.0001659
80 & 0.000671 & 0.00020130000000000001
90 & 0.000782 & 0.0002346
100 & 0.000935 & 0.0002805
110 & 0.00107 & 0.000321
120 & 0.00122 & 0.00036599999999999995
    """

    params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VT/pG * linfit_params[0]
    #print("TP Leck 5E-5 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $5 \cdot 10^{-5}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("TP_Leck_5e5.pdf")
    plt.figure()

    # # 7e-5 mbar
    pG = unp.uarray(7E-5, (7E-5)*0.3)
    t, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/TP_Leck_7e5_D2.csv", delimiter=",", skip_header=1, usecols=(2,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D2_err = D2*0.3
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung")

    #print("TP Leckrate 7e-5 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D2[index]} & {D2_err[index]}")
    #print(" ")
    """
    TP Leckrate 7e-5 mbar
t     p(t)     err
0 & 6.88e-05 & 2.0640000000000002e-05
10 & 0.000176 & 5.2799999999999996e-05
20 & 0.000285 & 8.549999999999999e-05
30 & 0.000384 & 0.0001152
40 & 0.000483 & 0.0001449
50 & 0.00061 & 0.00018299999999999998
60 & 0.000789 & 0.00023669999999999998
70 & 0.00101 & 0.000303
80 & 0.0012 & 0.00035999999999999997
90 & 0.00143 & 0.000429
100 & 0.00165 & 0.000495
110 & 0.00192 & 0.000576
120 & 0.00219 & 0.000657
    """

    params, cov = curve_fit(linfit, t[40:], D2[40:], sigma=D2_err[40:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit")
    S = VT/pG * linfit_params[0]
    #print("TP Leck 7E-5 mbar S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für $p_G$ = $7 \cdot 10^{-5}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("TP_Leck_7e5.pdf")
    plt.figure()

    ### Leitwerte
    # t, D1, D2 = usecols(2, 4, 5)

    # Fehler für Evakuierung:
    # 100 - 10^-8 mbar: 30%

    ### Gleichgewichtsdruck 2E-4 mbar
    pG = unp.uarray(2E-4, (2E-4)*0.3)
    t, D1, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/Leitparam_2e4_D1D2.csv", delimiter=",", skip_header=1, usecols=(2,4,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D1_err = D2*0.3
    D2_err = D2*0.3
    plt.errorbar(t, D1, yerr=D1_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung D1")
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="yellow", errorevery=(4), label="Messung D2")
    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für die Leitparameter bei $p_G$ = $2 \cdot 10^{-4}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("Leitparam_2e4.pdf")
    plt.figure()

    #print("Leitwert 2e-4 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D1[index]} & {D1_err[index]} & {D2[index]} & {D2_err[index]}")
    #print(" ")

    """
    Leitwert 2e-4 mbar
t     p(t)     err
0 & 0.0002 & 0.000429 & 0.00143 & 0.000429
10 & 0.00293 & 0.000849 & 0.00283 & 0.000849
20 & 0.0045 & 0.001338 & 0.00446 & 0.001338
30 & 0.00575 & 0.0018149999999999998 & 0.00605 & 0.0018149999999999998
40 & 0.00764 & 0.002346 & 0.00782 & 0.002346
50 & 0.00961 & 0.002979 & 0.00993 & 0.002979
60 & 0.0109 & 0.00339 & 0.0113 & 0.00339
70 & 0.012 & 0.00375 & 0.0125 & 0.00375
80 & 0.0134 & 0.00411 & 0.0137 & 0.00411
90 & 0.0147 & 0.0045 & 0.015 & 0.0045
100 & 0.0164 & 0.004889999999999999 & 0.0163 & 0.004889999999999999
110 & 0.0182 & 0.005399999999999999 & 0.018 & 0.005399999999999999
120 & 0.0199 & 0.0059099999999999995 & 0.0197 & 0.0059099999999999995
130 & 0.0212 & 0.006449999999999999 & 0.0215 & 0.006449999999999999
140 & 0.0226 & 0.00684 & 0.0228 & 0.00684
150 & 0.024 & 0.0072299999999999994 & 0.0241 & 0.0072299999999999994
160 & 0.0254 & 0.007619999999999999 & 0.0254 & 0.007619999999999999
170 & 0.0269 & 0.00804 & 0.0268 & 0.00804
180 & 0.0286 & 0.00852 & 0.0284 & 0.00852
190 & 0.0302 & 0.008969999999999999 & 0.0299 & 0.008969999999999999
200 & 0.032 & 0.00942 & 0.0314 & 0.00942
    """
    ### Gleichgewichtsdruck 5e-5 mbar
    t, D1, D2 = np.genfromtxt("Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/Leitparam_5e5_D1D2.csv", delimiter=",", skip_header=1, usecols=(2,4,5), dtype=None, encoding=None, unpack=True) # Unpack Data
    t = t - t[0]
    D1_err = D2*0.3
    D2_err = D2*0.3
    plt.errorbar(t, D1, yerr=D1_err, fmt=".", ecolor="lightblue", errorevery=(5), label="Messung D1")
    plt.errorbar(t, D2, yerr=D2_err, fmt=".", ecolor="yellow", errorevery=(4), label="Messung D2")

    #print("Leitwert 2e-4 mbar")
    #print("t     p(t)     err")
    #for index, time in enumerate(t):
    #    if time % 10 == 0:
    #        print(f"{time} & {D1[index]} & {D1_err[index]} & {D2[index]} & {D2_err[index]}")
    #print(" ")

    """
    Leitwert 5e-5 mbar
t     p(t)     err
0 & 4.95e-05 & 6.69e-05 & 0.000223 & 6.69e-05
10 & 0.000393 & 8.609999999999999e-05 & 0.000287 & 8.609999999999999e-05
20 & 0.000503 & 0.0001059 & 0.000353 & 0.0001059
30 & 0.00061 & 0.0001284 & 0.000428 & 0.0001284
40 & 0.00074 & 0.0001542 & 0.000514 & 0.0001542
50 & 0.000887 & 0.00018299999999999998 & 0.00061 & 0.00018299999999999998
60 & 0.00105 & 0.00022439999999999998 & 0.000748 & 0.00022439999999999998
70 & 0.00123 & 0.0002649 & 0.000883 & 0.0002649
80 & 0.00139 & 0.00030900000000000003 & 0.00103 & 0.00030900000000000003
90 & 0.00154 & 0.000348 & 0.00116 & 0.000348
100 & 0.00168 & 0.00039299999999999996 & 0.00131 & 0.00039299999999999996
110 & 0.00185 & 0.00044399999999999995 & 0.00148 & 0.00044399999999999995
120 & 0.00199 & 0.000501 & 0.00167 & 0.000501
130 & 0.00217 & 0.000552 & 0.00184 & 0.000552
140 & 0.00235 & 0.000606 & 0.00202 & 0.000606
150 & 0.00254 & 0.000651 & 0.00217 & 0.000651
160 & 0.00276 & 0.0007229999999999999 & 0.00241 & 0.0007229999999999999
170 & 0.00299 & 0.000795 & 0.00265 & 0.000795
180 & 0.00321 & 0.000849 & 0.00283 & 0.000849
190 & 0.00338 & 0.0009179999999999999 & 0.00306 & 0.0009179999999999999
200 & 0.00358 & 0.0009809999999999999 & 0.00327 & 0.0009809999999999999
    """

    # Saugvermögen für D1
    params, cov = curve_fit(linfit, t[20:], D1[20:], sigma=D1_err[20:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="red", label="Fit D1")
    S = VT/pG * linfit_params[0]
    #print("Leitwertmessung bei 7E-5 mbar D1 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

    # Saugvermögen für D2
    params, cov = curve_fit(linfit, t[20:], D2[20:], sigma=D2_err[20:], absolute_sigma=True)
    linfit_params = unp.uarray(params, np.sqrt(np.diag(cov)))
    x_fit = np.linspace(0, t[-1]+10)
    plt.errorbar(x_fit, linfit(x_fit, unp.nominal_values(linfit_params[0]), unp.nominal_values(linfit_params[1])), color="green", label="Fit D2")
    S = VT/pG * linfit_params[0]

    plt.xlabel("Zeit $t$ in $\si{\second}$")
    plt.ylabel("Druck $p$ in $\mathrm{mbar}$")
    plt.title("Leckratenmessung für die Leitparameter bei $p_G$ = $5 \cdot 10^{-5}$ $\mathrm{mbar}$")
    plt.legend()
    plt.savefig("Leitparam_5e5.pdf")
    #print("Leitwertmessung bei 7E-5 mbar D2 S = ", S, " l/s")
    #print("m = ", linfit_params[0], "   b = ", linfit_params[1])
    #print(" ")

'''
TP Fit 1 S =  7.6+/-0.8  l/s
m =  -0.2304202+/-0.0000004    b =  -2.030883+/-0.000005
 
TP Fit 2 S =  0.78+/-0.08  l/s
m =  -0.023777365+/-0.000000013    b =  -5.2276858+/-0.0000005
 
TP Fit 3 S =  0.192+/-0.019  l/s
m =  -0.0058313180+/-0.0000000027    b =  -6.12055098+/-0.00000026
 
TP Leck 1E-4 mbar S =  15+/-5  l/s
m =  (4.40+/-0.32)e-05    b =  -0.00097+/-0.00021
 
TP Leck 2E-4 mbar S =  23+/-8  l/s
m =  0.000139+/-0.000012    b =  -0.0012+/-0.0008
 
TP Leck 5E-5 mbar S =  7.5+/-2.4  l/s
m =  (1.14+/-0.08)e-05    b =  -0.00019+/-0.00006
 
TP Leck 7E-5 mbar S =  10.0+/-3.2  l/s
m =  (2.12+/-0.14)e-05    b =  -0.00045+/-0.00010
 
Leitwertmessung bei 7E-5 mbar D1 S =  2.7+/-0.8  l/s
m =  (1.61+/-0.06)e-05    b =  0.00011+/-0.00004
 
Leitwertmessung bei 7E-5 mbar D2 S =  2.4+/-0.8  l/s
m =  (1.44+/-0.06)e-05    b =  (-4+/-4)e-05
'''