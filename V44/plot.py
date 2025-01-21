import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import cmath
from tqdm import tqdm
import optuna
import uncertainties as uc
import uncertainties.unumpy as unp


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Be careful with rad/deg
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def mygenfromtxt(filename):
    header_lines = sum(1 for line in open(filename) if line.strip().startswith(';') or line.strip().startswith('_'))
    return np.genfromtxt(filename, skip_header=header_lines, unpack=True)

def gaussian(x, mean, sigma,a,b):
    return b+a*1.0 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

def gaussian_error(x, mean, sigma,a,b):
    return b+a*1.0 / (sigma * np.sqrt(2 * np.pi)) * unp.exp(-0.5 * ((x - mean) / sigma) ** 2)

def rel(x_true,x_meas):
    return abs(x_true - x_meas) / x_true *100

def calc_thickness(x,y,lamda):
    minima, _ = find_peaks(-y)
    minima=minima[4:11]
    minima_distances = [abs(x[minima[i]]-x[minima[i+1]]) for i in range(len(minima) - 1)]
    mean_kiessig_angle=np.mean(minima_distances)
    mean_kiessig_angle=uc.ufloat(np.deg2rad(mean_kiessig_angle),np.deg2rad(err_ofmean(minima_distances)))
 
    return lamda/(2*mean_kiessig_angle), minima

def err_ofmean(x):
    return np.std(x, ddof=1) / np.sqrt(len(x))

#First Plot 
# measured reflectivity, diffuse reflectivity, corrected reflectivity
ms_angle,ms_reflectivity=mygenfromtxt('data/Reflect1.UXD')
dif_angle,dif_reflectivity=mygenfromtxt('data/Reflect2.UXD')
ms_minus_dif=ms_reflectivity-dif_reflectivity



#create detector scan with gaussian fit
detector_scan_angle,detector_scan_reflectivity=mygenfromtxt(f'data/Detector1.UXD')

params_gaussian, pcov_gussian =curve_fit(gaussian,detector_scan_angle,detector_scan_reflectivity,p0=[0,0.05,1e5,0])
x=np.linspace(-0.5,0.5,1000)
y=gaussian_error(x,
           uc.ufloat(params_gaussian[0], np.sqrt(pcov_gussian[0,0])),
           uc.ufloat(params_gaussian[1], np.sqrt(pcov_gussian[1,1])),
           uc.ufloat(params_gaussian[2], np.sqrt(pcov_gussian[2,2])),
           uc.ufloat(params_gaussian[3], np.sqrt(pcov_gussian[3,3])))

max_intensity=np.max(y)
max_intensity_with_error=max_intensity
max_intensity=max_intensity.n
half_max = max_intensity / 2
indices=np.where(y >= half_max)

fwhm=abs(x[indices[0][-1]]-x[indices[0][0]])

plt.figure()
plt.plot(ms_angle,ms_reflectivity/(5*max_intensity) ,'.',markersize=1,label='Measured Reflectivity')
plt.plot(ms_angle,dif_reflectivity/(5*max_intensity),'.',markersize=1,label='Diffuse Reflectivity')
plt.plot(ms_angle,ms_minus_dif/(5*max_intensity)    ,'.',markersize=1,label='Measured-Diffuse')
plt.xlabel(r'$\alpha_\text{i}\,/\,°$')
plt.ylabel('Reflectivity')
plt.yscale('log')
plt.legend()
plt.grid()
plt.savefig('build/basic_reflectivity_curves.pdf')

plt.figure()
plt.plot(detector_scan_angle,detector_scan_reflectivity,'x',label='Detector Scan')
plt.plot(x,unp.nominal_values(y),label='fit')
plt.axvline(x[indices[0][0]],color='black', linestyle='--',label='FWHM')
plt.axvline(x[indices[0][-1]],color='black', linestyle='--')
plt.xlabel(r'$\alpha\,/\,°$')
plt.ylabel('I')
plt.legend()
plt.grid()
plt.savefig(f'build/detector_scan_1.pdf')


#Second Plot
# Diffuse-Measured, ideal fresnel reflectivity , corrected by geometry factor

#Ideal Fresnel Reflectivity
# ------------------------------------------
def R_fresnel_approx(ac,ai):
    return (ac/(2*ai+1e-10))**4

# ac=0.233
# ac_index=np.abs(ms_angle - ac).argmin()
# ac_index=0
# fresnel_angle=ms_angle[ac_index:]
# ideal_fresnel_plustotalreflexion=np.concatenate([np.ones_like(ms_angle[:ac_index]),ideal_fresnel])

def R_fresnel(ai):
    wavelength=1.541e-10
    k = 2 * np.pi / wavelength
    n1=1 
    delta2=7.6e-06
    beta2=(1.541e-10/(4*np.pi))*141*1/(100)
    n2=1-delta2-1j*beta2    
    kz1 = k*np.sqrt(n1**2-np.cos(np.deg2rad(ai))**2)
    kz2 = k*np.sqrt(n2**2-np.cos(np.deg2rad(ai))**2)

    r = ((kz1-kz2)/(kz1+kz2))
    return abs(r)**2

ideal_fresnel=R_fresnel(ms_angle)

# ------------------------------------------

# Correction by geometry factor
# ------------------------------------------

#calculate beam width
zscan1_angle,zscan1_reflectivity=mygenfromtxt(f'data/Zscan1.UXD')

last_max_index = np.where(np.logical_and((np.max(zscan1_reflectivity)-10000 <=zscan1_reflectivity),(zscan1_reflectivity <= np.max(zscan1_reflectivity))))[0][-1]
first_min_index = np.where(np.logical_and((np.min(zscan1_reflectivity)+10000 >=zscan1_reflectivity),(zscan1_reflectivity >= np.min(zscan1_reflectivity))))[0][0]

beam_width=abs(zscan1_angle[first_min_index]-zscan1_angle[last_max_index])

plt.figure()
plt.plot(zscan1_angle,zscan1_reflectivity,'x',label='Z Scan')
plt.axvline(zscan1_angle[last_max_index],label='lastmax',color='black', linestyle='--')
plt.axvline(zscan1_angle[first_min_index],label='firstmin',color='black', linestyle='--')
plt.xlabel(r'$z$\,/\,mm')
plt.ylabel('$I$')
plt.legend()
plt.grid()
plt.savefig(f'build/z_scan_1.pdf')

D=2e-2
d0=beam_width*1e-3

def geomatry_factor(ai,ag):
    geometry_factor=[]
    for i in ai:
        if 0<i<ag:
            geometry_factor.append((D*np.sin(np.deg2rad(i)))/(d0))
        else:
            geometry_factor.append(1)
    return np.array(geometry_factor)
theoretical_geometry_angle=np.rad2deg(np.arcsin(d0/D))

rockingscan_angle,rockingscan_reflectivity=mygenfromtxt(f'data/Rocking1_2.UXD')
maxindex=np.where(rockingscan_reflectivity==np.max(rockingscan_reflectivity))[0][0]

def linear(x,a,b):
    return a*x+b
def xzero(a,b):
    return -b/a

params_left,pcov_left=curve_fit(linear,rockingscan_angle[0:maxindex],rockingscan_reflectivity[0:maxindex])
params_right,pcov_right=curve_fit(linear,rockingscan_angle[maxindex:],rockingscan_reflectivity[maxindex:])

x_left=np.linspace(rockingscan_angle[0],rockingscan_angle[maxindex],1000)
x_right=np.linspace(rockingscan_angle[maxindex],rockingscan_angle[-1],1000)
xzero_left=xzero(params_left[0],params_left[1])
xzero_right=xzero(params_right[0],params_right[1])
measured_geometry_angle=abs(xzero_left-xzero_right)/2

plt.figure()
plt.plot(rockingscan_angle,rockingscan_reflectivity,'x',label='Rocking Scan')
plt.plot(x_left,linear(x_left,params_left[0],params_left[1]))
plt.plot(x_right,linear(x_right,params_right[0],params_right[1]))
plt.xlabel(r'$y$\,/\,mm')
plt.ylabel('$I$')
plt.legend()
plt.grid()
plt.savefig(f'build/rockingscan.pdf')


geometry_corrected_relativity=ms_minus_dif*(1/geomatry_factor(ms_angle,measured_geometry_angle))
# ------------------------------------------

plt.figure()
plt.plot(ms_angle, ms_minus_dif/(5*max_intensity),'.',markersize=1, label='Measured-Diffuse')
plt.plot(ms_angle, ideal_fresnel,label='Ideal Fresnel Reflectivity')
# plt.plot(fresnel_angle, ideal_fresnel,label='Ideal Fresnel Reflectivity')
plt.plot(ms_angle,geometry_corrected_relativity/(5*max_intensity),'.',markersize=1,label='Corrected by geometry factor')
plt.xlabel(r'$\alpha_\text{i}\,/\,°$')
plt.ylabel('Reflectivity')
plt.yscale('log')
plt.legend()
plt.grid()
plt.savefig('build/advanced_curves.pdf')


# Third Plot
# corrected by geometry factor, by parrat algorithm 
# from this inferre -> dispersion, thikness-> critical angle

def distance(calced_dist,dist_comp):  
    result=np.sum(np.abs(np.log(1.0/(calced_dist+1e-6))-np.log(1.0/(dist_comp+1e-6))))*(1/len(calced_dist))  
    return result

def optuna_objective(trial,dist_comp,parrat_angle):
    layer_thickness = trial.suggest_float('layer_thickness', 8.5e-8,8.7e-8)
    delta1  = trial.suggest_float('delta1', 3e-6,8e-6 )
    delta2  = trial.suggest_float('delta2', 8e-7,3e-6 )
    beta1   = trial.suggest_float('beta1', 6e-8,3e-7  )
    beta2   = trial.suggest_float('beta2', 1e-9,9e-9  )
    sigma1  = trial.suggest_float('sigma1', 1e-11,1e-9)
    sigma2  = trial.suggest_float('sigma2', 1e-11,1e-9)
    return distance(gen_parrat_curve(parrat_angle,1.54e-10,[layer_thickness,0],[delta1,delta2,0],[beta1,beta2,0],[sigma1,sigma2]),dist_comp)

def parameter_search(n_trials,dist_comp,parrat_angle):
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(study_name='parratt',sampler=(optuna.samplers.TPESampler(seed=1)))
    study.optimize(lambda trial: optuna_objective(trial,dist_comp,parrat_angle), n_trials=n_trials, catch=True,show_progress_bar=True)
    return study.best_params


def parratt_reflectivity(angle, wavelength, layer_thickness, num_layers,delta,beta,sigma):
    X = np.zeros(num_layers+1,dtype=complex)
    r = np.zeros(num_layers+1, dtype=complex)
    k = 2 * np.pi / wavelength
    kz= np.zeros(num_layers+1,dtype=complex)
    n = np.zeros(num_layers+1,dtype=complex)
    n[0]=1-delta[0]-1j*beta[0]
    kz[0] = k*np.sqrt(n[0]**2-np.cos(np.deg2rad(angle))**2)

    for i in range(num_layers):
        n[i+1]=1-delta[i+1]-1j*beta[i+1]     
        kz[i+1] = k*np.sqrt(n[i+1]**2-np.cos(np.deg2rad(angle))**2)
        r[i] = ((kz[i+1]-kz[i])/(kz[i+1]+kz[i])) *  np.exp(-2*kz[i]*kz[i+1]*sigma[i]**2)
        X[i+1]=np.exp(-2j*kz[i+1]*layer_thickness[i])*((r[i] + X[i])/(1+r[i]*X[i] ))       
    return abs(X[-1])**2

def gen_parrat_curve(parrat_angle, wavelength, layer_thickness, delta, beta, sigma):
    parrat_reflectivity=np.empty(len(parrat_angle))
    for index,angle in enumerate(parrat_angle):
        parrat_reflectivity[index] = parratt_reflectivity(angle, wavelength, layer_thickness, len(layer_thickness),delta,beta,sigma)
    return parrat_reflectivity
    
#parameters from the altprotocoll
wavelength = 1.54e-10  # Wavelength in Angstroms
layer_thickness = [8.6e-8,0]  # Layer thicknesses in Angstroms
delta=[6e-6, 0.6e-6, 0]
beta=[delta[0]/200,delta[1]/40,0]
sigma=[6.45*10**(-10),5.5*10**(-10)]#roughness

#altprotcoll result
parrat_angle=np.linspace(0, 2.5, 2000)
parratt_refl = gen_parrat_curve(parrat_angle, wavelength, layer_thickness, delta, beta, sigma)


#perform parametersearch using optuna
parrat_angle_optuna = np.linspace(0, 2.5, len(geometry_corrected_relativity))
optuna_params=parameter_search(300,geometry_corrected_relativity[1:]/(5*max_intensity),parrat_angle_optuna[1:])
parratt_refl_optuna = gen_parrat_curve(parrat_angle, 1.54e-10,[optuna_params['layer_thickness'],0],
                                                           [optuna_params['delta1'],optuna_params['delta2'],0],
                                                           [optuna_params['beta1'],optuna_params['beta2'],0],
                                                           [optuna_params['sigma1'],optuna_params['sigma2']])

#calc layer thickness
thickness_simple_background_correction,minima_simple=calc_thickness(ms_angle,ms_minus_dif,1.54e-10)
thickness_geometry_correction,minima_geometry=calc_thickness(ms_angle,geometry_corrected_relativity,1.54e-10)
thickness_parrat,minima_parrat=calc_thickness(parrat_angle,parratt_refl_optuna,1.54e-10)

plt.figure()
plt.plot(ms_angle[1:], geometry_corrected_relativity[1:]/(5*max_intensity),label='Corrected by geometry factor')
#for i in minima_geometry:
#    plt.axvline(ms_angle[i],color='black', linestyle='--')
#plt.axvline(ms_angle[minima_geometry[0]],color='black', linestyle='--',label='minima')
plt.plot(parrat_angle, parratt_refl_optuna ,label='Parratt')
plt.xlabel(r'$\alpha_\text{i}\,/\,°$')
plt.ylabel('Reflectivity')
plt.yscale('log')
plt.legend()
plt.grid()
plt.savefig('build/parrat_comparison.pdf')


#write results to txt file
with open('build/calculations.txt', 'w') as file:
    file.write(f'Detector Scan:\n')
    file.write(f'Parameters of gaussian fit:\n')
    file.write(f'mu=   {uc.ufloat(params_gaussian[0], np.sqrt(pcov_gussian[0,0]))}\n')
    file.write(f'sigma={uc.ufloat(params_gaussian[1], np.sqrt(pcov_gussian[1,1]))}\n')
    file.write(f'a=    {uc.ufloat(params_gaussian[2], np.sqrt(pcov_gussian[2,2]))}\n')
    file.write(f'b=    {uc.ufloat(params_gaussian[3], np.sqrt(pcov_gussian[3,3]))}\n')
     
    file.write(f'Full width half height: {fwhm}\n')
    file.write(f'Maximum Intensity: {max_intensity_with_error}\n')
    file.write(f'\n')
    file.write(f'Z-Scan\n')
    file.write(f'Beam width: {beam_width}\n')
    file.write(f'\n')
    file.write(f'Rockingscan\n')
    file.write(f'Theoretical geometry angle: {theoretical_geometry_angle}\n')
    file.write(f'Measured geometry angle: {measured_geometry_angle}\n')
    file.write(f'\n')
    file.write(f'Layer Thickness with simple background correction: {thickness_simple_background_correction}\n')
    file.write(f'Layer Thickness with geometry correction: {thickness_geometry_correction}\n')
    file.write(f'Layer Thickness from parratt: {thickness_parrat}\n')
    file.write(f'Parameters from parrat and optuna/parametersearch: \n{optuna_params}\n')
    ac_si=np.rad2deg(np.sqrt(2*optuna_params["delta1"]))
    ac_pol=np.rad2deg(np.sqrt(2*optuna_params["delta2"]))
    file.write(f'Critical angle (Si): {ac_si}\n')
    file.write(f'Critical angle (Pol): {ac_pol}\n')
    file.write(f'\n')
    file.write(f'Discussion\n')
    file.write(f'rel ac lit si: {rel(0.223,ac_si)}\n')
    file.write(f'rel ac lit pol: {rel(0.153,ac_pol)}\n')
    file.write(f'rel layer thickness kiessig parratt: {rel(thickness_geometry_correction,optuna_params["layer_thickness"])}\n')
    file.write(f'rel dispersion parratt lit si: {rel(7.6e-6,optuna_params["delta1"])}\n')
    file.write(f'rel dispersion parratt lit pol: {rel(3.5e-6,optuna_params["delta2"])}\n')
    file.write(f'rel geometry angle: {rel(theoretical_geometry_angle,measured_geometry_angle)}\n')
    file.write(f'\n')
    file.write(f'\n')
    file.write(f'\n')
    file.write(f'\n')


# preview plots
files=['rocking_scan_1_1',
       'rocking_scan_1_2',
       'rocking_scan_1',
       'rocking_scan_2_better',
       'rocking_scan_2',
       'rocking_scan_3',
       'z_scan_1',
       'z_scan_2',
       'z_scan_3',
       'x_scan_1',
       'detector_scan_1']
for i in files:
    angle,reflectivity=mygenfromtxt(f'datastefan/{i}.UXD')

    plt.figure()
    plt.plot(angle,reflectivity,'x',label='Rocking Scan')
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Reflectivity')
    plt.legend()
    plt.grid()
    plt.savefig(f'build/preview_{i}.pdf')






















#Comments contain code for a gridsearch and a curvefit

# def generate_combinations(arrays, current_combination, index):
#     if index == len(arrays):
#         yield tuple(current_combination)
#     else:
#         for item in arrays[index]:
#             current_combination.append(item)
#             yield from generate_combinations(arrays, current_combination, index + 1)
#             current_combination.pop()

# def grid_search(parameter_combinations,dist_comp,objective):     
#        best_params = None
#        best_objective = float('inf') 
#        for params in tqdm(parameter_combinations):     
#            objective_value = distance(objective(len(dist_comp),1.54e-10,[params[0],0],[params[1],params[2],0],[params[1]/200,params[2]/40,0],[params[3],params[4]],)[1],dist_comp)
#            if objective_value < best_objective:
#                best_objective = objective_value
#                best_params = params   
#        return best_params                    
                      
# def gridsearch_plus_parameterhandling(precision,i):
#     best_params=[] 
#     #getting prameterspaces
#     searchspaces={"1": {"layer_thickness": [8.4e-8,8.8e-8], 
#                         "delta1":          [4e-6,8e-6],   
#                         "delta2":          [0.1e-6,1e-6], 
#                         # "beta1":         [1e-8,5e-7],   
#                         # "beta2":         [1e-10,1e-9],   
#                         "sigma1":          [5e-10,8e-10],  
#                         "sigma2":          [4e-10,7e-10]}}
#     for key in searchspaces[i].keys():
#         searchspaces[i][key]=np.linspace(searchspaces[i][key][0],searchspaces[i][key][1],precision)
 
#     array_list = [value for key,value in searchspaces[i].items()]
#     parameter_combinations = list(generate_combinations(array_list, [], 0))
   
#     best_params=grid_search(parameter_combinations,ms_minus_dif,gen_parrat_curve)          
#     return best_params



# precision=8
# best_params=gridsearch_plus_parameterhandling(precision,'1')
# parrat_angle_search,parratt_refl_search = gen_parrat_curve(2000, 1.54e-10,[best_params[0],0],[best_params[1],best_params[2],0],[best_params[1]/40,best_params[2]/40,0],[best_params[3],best_params[4]])
# print(f'precision {precision}\n{best_params}')

# def gen_parrat_curvefit(parrat_angle, layer_thickness, delta1, delta2, beta1, beta2, sigma1, sigma2):
#     wavelength=1.54e-10
#     parrat_reflectivity=np.empty(len(parrat_angle))
#     for index,angle in enumerate(parrat_angle):
#         parrat_reflectivity[index] = parratt_reflectivity(angle, wavelength, [layer_thickness,0], 2 ,[delta1,delta2,0],[beta1,beta2,0],[sigma1,sigma2])
#     return parrat_reflectivity

# parrat_angle_curvefit = np.linspace(0, 2, len(ms_reflectivity))
# print(len(parrat_angle_curvefit))
# params_parrat_curvefit,_=curve_fit(gen_parrat_curvefit,parrat_angle_curvefit[50:200], ms_reflectivity[50:200]/(5*max_intensity),p0=[8.6e-8,6e-6,0.6e-6,6e-6/200,0.6e-6/40,6.45*10**(-10),5.5*10**(-10)],maxfev=2000)

# parrat_angle_curvefit, parratt_refl_curvefit = gen_parrat_curve(2000, 1.54e-10,[params_parrat_curvefit[0],0],
#                                                            [params_parrat_curvefit[1],params_parrat_curvefit[2],0],
#                                                            [params_parrat_curvefit[3],params_parrat_curvefit[4],0],
#                                                            [params_parrat_curvefit[5],params_parrat_curvefit[6]],)
# print(params_parrat_curvefit)













# def parratt_reflectivity(angle, wavelength, layer_thickness, num_layers,delta,beta,sigma):
#     r = np.zeros(num_layers+1, dtype=complex)
#     R = np.zeros(num_layers+1,dtype=complex)
#     T = np.ones(num_layers+1,dtype=complex)
#     k = 2 * np.pi / wavelength
#     kz = np.zeros(num_layers+1,dtype=complex)
#     n = np.zeros(num_layers+1,dtype=complex)

#     for i in range(num_layers):
#         #komplexe zahlen unter wurzel prüfen
#         # vorzeichen bei n überprüfen
#         n[i]=1-delta[i]+1j*beta[i]
#         n[i+1]=1-delta[i+1]+1j*beta[i+1]

#         kz[i] = k*np.sqrt(n[i]**2-np.cos(np.deg2rad(angle))**2)
#         kz[i+1] = k*np.sqrt(n[i+1]**2-np.cos(np.deg2rad(angle))**2)
          
#         r[i] = ((kz[i+1]-kz[i])/(kz[i+1]+kz[i])) *  np.exp(-2*kz[i]*kz[i+1]*sigma[i]**2)

#         R[i+1]=(1.0/(1-r[i])) * ((T[i])*r[i]*np.exp(-1j*(kz[i]+kz[i+1])*layer_thickness[i])+R[i]*np.exp(-1j*(kz[i]-kz[i+1])*layer_thickness[i]))
#         T[i+1]=(1.0/(1-r[i])) * ((T[i])*     np.exp(1j*(kz[i]-kz[i+1])*layer_thickness[i]) +R[i]*r[i]*np.exp(1j*(kz[i]+kz[i+1])*layer_thickness[i]))
#     return abs(R[-1])**2
