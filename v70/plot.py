import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

#### Drehschieberpumpe

# t und D0 = usecols(2, 3)

# Fehler für Evakuierung
# 1000 - 10 mbar: 0.3%
# 10 - 2*10^-3:    10%

# Fehler für Leckratenmessung
# 0.5 mbar:        10%
# Sonst:          0.3%

###################

#### Turbomolekularpumpe

# t und D2 = usecols(2, 5)

# Fehler für Evakuierung
# 100 - 10^-8 mbar: 30%

# Fehler für Leckratenmessung
# 100 - 10^-8 mbar: 30%

### Leitwerte
# t, D1, D2 = usecols(2, 4, 5)

# Fehler für Evakuierung:
# 100 - 10^-8 mbar: 30%

