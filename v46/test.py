import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl

me = constants.electron_mass
m_lit = 1.9*me
m_ex = np.array([2.9e-31, 3.9e-31, 1.6e-31, 2.2e-31])
diff = m_ex - m_lit
print(diff/m_lit)