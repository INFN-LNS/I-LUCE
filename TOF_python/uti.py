import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import os


#read file with ToF data from scope
def read_tof(f_name:str):

    #chek if file_path exists:
    if os.path.exists(f_name) == False:
        print('File does not exist')
    else:
        data=pd.read_csv(f_name,header=19,sep=',',engine='python')
        data=np.array(data,dtype=float)

    time,ch1,ch2,ch3,ch4 = np.split(data,5,axis=1)

    return time,ch1,ch2,ch3,ch4


#functio to read SRIM output file on stopping power and give two np arrays: incident energy in keV and total stopping power for that energy in keV/um

def read_SRIM_output(f_name:str):

    #chek if file_path exists:
    if os.path.exists(f_name) == False:
        print('File does not exist')

    else:
            
        d = pd.read_csv(f_name, header=20,skipfooter=13,sep=r'\s+',engine='python')
        d = np.array(d)
    
        n_rows, n_columns = d.shape
        
        E = []
        total_stopping_power = []
        Range = []
        
        for i in range(0,n_rows):
            d[i][0] = float(d[i][0].replace(',', '.'))
            d[i][2] = float(d[i][2].replace(',', '.'))
            d[i][3] = float(d[i][3].replace(',', '.'))
            d[i][4] = float(d[i][4].replace(',', '.'))
            if d[i][1] == 'MeV':
                d[i][0] = d[i][0]*1000
            if d[i][5] == 'A':
                d[i][4] = d[i][4]*0.0001
            
            E.append(d[i][0])
            total_stopping_power.append(d[i][2]+d[i][3])
            Range.append(d[i][4])

    E=np.array(E)
    total_stopping_power=np.array(total_stopping_power)
    Range=np.array(Range)

    return E,total_stopping_power,Range



#function to calculate energy loss of particle with energy E0 [keV] inside a target of thickness d [um]. Particle tipe (proton, alpha ecc.) and target properties are given with SRIM_ouput file on stopping power  

def e_loss(f_name:str,E0,d,n_iter):

    loss = 0
    #chek if file_path exists:
    if os.path.exists(f_name) == False:
        print('File does not exist')

    else: 

        E,sp,_ = read_SRIM_output(f_name)
        interp_dEdX = interp1d(E,sp, kind ='linear')

        dx = d/n_iter
        e_in = np.copy(E0)
        #loss = 0
      
        
        for i in range (0, n_iter):
            loss_step = interp_dEdX(e_in)*dx
            e_in -= loss_step
            if (e_in <= 2).any():
                break
            else:
                loss += loss_step
           

    return loss

#function to get the max range of a particle with energy E_in [keV] in a target. Partcle and target properties are given with the SRIM outpu file
def particle_range(f_name:str,E_in):

    #chek if file_path exists:
    if os.path.exists(f_name) == False:
        print('File does not exist')

    else:

        E,_,r = read_SRIM_output(f_name)
        interp_range = interp1d(E,r, kind ='linear')


    return interp_range(E_in)



#function to compute the shifted Maxwell-Boltzman equation signal for Time Of Flight: parametri(tempo,distanza,massa ione, energia termica Kb*T)
def SMB_in_time(t,L,m,kT,vs):
    
    vp = L/t
    
    normalization_constant = 2*pow(np.pi,3)*m/kT
    normalization_constant = pow(normalization_constant,1/2)
    
    beta = m/(2*kT)

    distribution = pow(beta/np.pi,3/2)
    distribution = distribution*pow(L,4)/pow(t,5)
    distribution = distribution*np.exp(-beta*pow(vp-vs,2))

    return normalization_constant*distribution

    
#class to construct detector:

class Detector:
    def __init__(self,material,area,thickness):
        self.material = material
        self.area = area
        self.thickness = thickness
        if self.material == 'Diamond':
            self.e_h = 13 #eV
            self.density = 3.51 # g/cm^-3
        if self.material == 'SiC':
            self.e_h =  7.6#eV
            self.density = 3.21 # g/cm^-3









