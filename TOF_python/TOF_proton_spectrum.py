#import ROOT
import uti
#from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np



#Physycs constant:

mp = 1.627e-27 #proton mass [kg]
q = 1.602e-19 #elementary charge [C]
c =  3e8 #light speed

##########################################


#Define all the experimental parameters needed for the ToF signal processing:

diamond1 = uti.Detector('Diamond',6.3,35) #detector definition: material('Diamond' or 'SiC'), Sensitive area [mm^2], thickness [um]
L = 1.07 #distance between target and detector [m]
Al_filter_thickness = 2.7 #Aluminum filter thikness in um
R = 50
omega = (diamond1.area*1e-6)/pow(L,2)

##########################################



#Initialize file for data analysis:

shot_number = 60715
file_name = f'/Users/mcsalmone/Desktop/PhD_Catania/activities/TOF/scope1/02_02_2024/{shot_number}.csv'
diamonds_scope_1 = uti.read_tof(file_name)

######################################################


time1 = diamonds_scope_1[0]
d1 = diamonds_scope_1[1]



# Graphic: draw the plot of the oscilloscope original ToF signal
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(time1, d1, label='Original ToF signal d1', marker='.', color='black')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Signal [V]')
ax.legend()
ax.legend()


cursor = Cursor(ax, horizOn=True, vertOn=True, useblit=True, color='red', linewidth=1)


#Photopeak selection from oscilloscope data:

selected_photopeak = None  # Stores the selected x-coordinate
print('Please select the photopeak time to translate the temporal data...(then close the window)')
def photopeak_selection(event):
    global selected_photopeak  # Ensure we modify the global variable
    if event.inaxes is not ax:
        return  # Ignore clicks outside the plot area
    
    selected_photopeak = event.xdata  # Store the selected x-coordinate
    print(f"Selected point: {selected_photopeak}")


fig.canvas.mpl_connect('button_press_event', photopeak_selection)
plt.show()


tof = time1 - selected_photopeak #move time signal so photopeak ia at time 0


# Graphic: draw the plot of the oscilloscope traslated ToF signal
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(tof, d1, label='ToF signal d1', marker='.', color='black')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Signal [V]')
ax.legend()
plt.title('ToF')

cursor = Cursor(ax, horizOn=True, vertOn=True, useblit=True, color='red', linewidth=1)



#Selection of the point for the signal of the fast protons:

t_min_proton = None # Stores the selected x-coordinate
print('Please select the fast proton time...(then close the window)')
def fast_proton_time_selection(event):
    global t_min_proton  # Ensure we modify the global variable
    if event.inaxes is not ax:
        return  # Ignore clicks outside the plot area
    
    t_min_proton = event.xdata  # Store the selected x-coordinate
    print(f"Selected point: {t_min_proton}")


fig.canvas.mpl_connect('button_press_event', fast_proton_time_selection)
plt.show()



v_proton_max = L/t_min_proton #Compute the maximum velocity for the proton an

e_max = 0.5*mp*pow(v_proton_max,2) #Compute the maximum proton energy
print('Maximum proton energy:',e_max*(1/q)*1e-6, 'MeV')

dt_only_proton = t_min_proton*(pow(2,0.5)-1)
t_max_only_proton = t_min_proton + dt_only_proton #Compute the maximum ToF signal in which there are only protons



#Now compute the energy spectrum in the ToF signal in which I am sure there are only protons:
spectrum_true = []
energy_proton = []
for s,t in zip (d1,tof):
    if (t >= t_min_proton) and (t <= t_max_only_proton):
        v_p = L/t
        e_p = 0.5*mp*pow(v_p,2)*(1/q)*1e-3
        energy_proton.append(e_p)
        e_loss_in_Al = uti.e_loss('SRIM_Outputs/Hydrogen in Aluminum.txt',e_p,Al_filter_thickness,1000)
        e_after_Al = e_p - e_loss_in_Al
        e_loss_in_Diamond = uti.e_loss('SRIM_Outputs/Hydrogen in Carbon.txt',e_after_Al,diamond1.thickness,1000)
        spectrum_true.append(0.5*diamond1.e_h*1e-3*s*t/(q*R*pow(e_loss_in_Diamond,2)*omega))





#Draw two plots: 1) TOF oscilloscop signal 2) differential proton spectrum
fig, ax = plt.subplots(2, figsize=(10, 7))

ax[0].plot(tof, d1, label='ToF signal d1', marker='.', color='black')
ax[0].axvline(x=t_min_proton, color='r', linestyle='--', label=f't_p_min = {t_min_proton}')
ax[0].axvline(x=t_max_only_proton, color='b', linestyle='--', label=f't_max_only_p = {t_max_only_proton}')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Signal [V]')
ax[0].legend()

ax[1].plot(energy_proton, spectrum_true, label='differential proton spectrum', marker='.', color='black')
ax[1].set_xlabel('E_p[keV]')
ax[1].set_ylabel(r'$\frac{d^{2}N}{dEd\Omega} [keV^{-1}sr^{-1}]$')
ax[1].set_yscale('log')
ax[1].legend()


#data = np.column_stack((energy_proton, spectrum_true))
#np.savetxt(f"{shot_number}.txt", data, fmt="%d", delimiter="\t", header="E_p\tSpectrum", comments="")

plt.show()











