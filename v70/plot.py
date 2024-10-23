import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl


# Welche Größen für die Auswertung?

#### Drehschieberpumpe

# t und D0 = usecols(2, 3)

# Fehler für Evakuierung
# 1000 - 10 mbar: 0.3%
# 10 - 2*10^-3:    10%

# Fehler für Belüftung
# 0.5 mbar:        10%
# Sonst:          0.3%

###################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
#
t, D0, D1, D2 = np.genfromtxt('Gruppe_Wieberg_Schink/FP_V70_Celina_Aaron/DP_Evak_D0.csv', delimiter=',', skip_header=1, usecols=(2, 3, 4, 5), dtype=None, encoding=None, unpack=True)

plt.plot(t, D0, ".", label="D0")
plt.plot(t, D1, ".", label="D1")
plt.plot(t, D2, ".", label="D2")
plt.yscale("log")
plt.legend()
plt.show()