import matplotlib.pyplot as plt         # Für Plots
import numpy as np                      # Fürs Rechnen
from uncertainties import ufloat        # zahl = ufloat(nominal_value, std_devs)
import uncertainties.unumpy as unp      # wie ufloat nur mit arrays
from scipy.optimize import curve_fit    # params, cov = curve_fit(fitfunktion, x-wert, y-wert, cov=True)
import scipy.constants as constants         # z.B. h = constants.h für planckzahl
import numpy.fft as fft
from skimage.feature import peak_local_max
from skimage.filters import gaussian

def sort_cyclic(points): # funktion um die peaks zyklisch zu sortieren

    pts = np.asarray(points) # makes sure dass wir ein array haben

    cx, cy = pts.mean(axis=0) # berechnet centroid

    angles = np.arctan2(pts[:,1] - cy, pts[:,0] - cx) # berechnet winkel 

    order = np.argsort(angles) # sortiert nach winkel

    return pts[order], order


px_imagesize = 128  # px
nm_imagesize = 2.31 # nm
scaling = nm_imagesize/px_imagesize # nm/px

show = False # parameter die bstrimmen ob wir bilder zeigen und oder speichern
save = False

def lattice_constant(number):
    data = gaussian(np.genfromtxt(f"data/{number}.csv", delimiter=";")) 

    plt.imshow(data, origin="lower", cmap="BuPu", extent=[0,nm_imagesize, 0,nm_imagesize]) # plot data
    plt.xlabel(r"$x$ in \si{\nano\meter}")
    plt.ylabel(r"$y$ in \si{\nano\meter}")
    plt.colorbar(label=r"Intensität in arbiträren Einheiten")
    plt.tight_layout()

    if show == True:
        plt.show()
    if save == True:
        plt.savefig(f"plots/{number}_meas.pdf")

    plt.figure()

    f = fft.ifftshift(data) # fouriertransformation in diesem abschnitt
    f = fft.fft2(f)
    f = fft.fftshift(f)
    f = np.abs(f)**2 # und power spectrum darstellen

    freqs = np.fft.fftfreq(px_imagesize, d=scaling) # so skalieren dass wir im richtigen reziproken raum (1/nm) sind
    freqs = np.fft.fftshift(freqs)
    kx = freqs
    ky = freqs

    peaks = peak_local_max(f, min_distance=5, num_peaks=6) # peaks finden und sortieren
    sorted_peaks, order = sort_cyclic(peaks)

    kx_axis = np.fft.fftshift(np.fft.fftfreq(px_imagesize, d=scaling))
    ky_axis = kx_axis.copy()

    kx_peaks = kx_axis[sorted_peaks[:,1]]
    ky_peaks = ky_axis[sorted_peaks[:,0]]
    

    plt.imshow(f, origin="lower", extent=[kx[0], kx[-1], ky[0], ky[-1]])
    plt.xlabel(r"$k_x$ in \si{\per\nano\meter}")
    plt.ylabel(r"$k_y$ in \si{\per\nano\meter}")
    plt.colorbar(label=r"Intensität in arbiträren Einheiten")
    plt.tight_layout()
    plt.colorbar()

    col = ["red", "blue", "green", "white", "orange", "purple"] # ignorieren
    #for i in range(6):
    #    plt.scatter(kx_peaks[i], ky_peaks[i], facecolors="none", edgecolors=col[i], s=100, label=f"{i}")
    #plt.legend()

    if show == True:
        plt.show()
    if save == True:
        plt.savefig(f"plots/{number}_ft.pdf")

    k_distances = np.sqrt((kx_peaks - np.roll(kx_peaks, -1))**2 + (ky_peaks - np.roll(ky_peaks, -1))**2) # berechne abstände zwischen den peaks
    k_distances = ufloat(np.mean(k_distances), np.std(k_distances)/np.sqrt(6)) 
    r_distances = 2/(np.sqrt(3)*k_distances) # berechne länge des gittervektors vom reziproken gittervektor
    print(f"Distances {number} = {r_distances} nm")

    plt.figure()

    return r_distances

save = True
r1 = lattice_constant(1)
lattice_constant(0)
r2 = lattice_constant(2)
r3 = lattice_constant(3)
r4 = lattice_constant(4)
r = np.array((unp.nominal_values(r1), unp.nominal_values(r2), unp.nominal_values(r3), unp.nominal_values(r4)))
r_avg = ufloat(np.mean(unp.nominal_values(r)), np.std(unp.nominal_values(r))/np.sqrt(4))
print("Average distance = ", r_avg, " nm")
print("Deviation from literature = ", 100*np.abs(r_avg - 0.246)/0.246, " %")



heights = np.array((1.64, 1.63, 1.11, 1.76, 1.35, 1.58))

height = ufloat(np.mean(heights), np.std(heights)/np.sqrt(len(heights)))

print(f"Height = {height} nm")