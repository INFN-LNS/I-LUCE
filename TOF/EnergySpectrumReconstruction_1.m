%% Experimental data 
% c'è un filtro e i protoni dopo il filtro si fermano dentro al rivelatore

% Rivelatore: SiC (densità=3.21 g/cm3)
% Spessore Rivelatore:27 um
% Area Rivelatore= 1 mm2
% Distanza target-rivelatore= 113.5 cm

% Spessore filtro alluminio=6.5 um Al
% Attenuatore=20 DB

data = textread('/Users/Ferdi/OneDrive/Documenti/MATLAB/DatiProva Spettro/ConFiltroCheSiFermano/C1shot100000.txt', '', 'delimiter', '.', ... 
                'headerlines',6,'emptyvalue', NaN);

% Commenta SiC o Diamante e controlla le aree S_d
           S='SiC';
           %S='Diamond'; 
           dB=10; %%-> Attenuation 20 DB->10; 10 DB->3; 6 DB->2; 13 DB-> 4.45
  
if  strcmp(S,'SiC')
   
figure;
plot(data(:,1)*10^9,smooth(data(:,2)*dB,20));
xlim([8000 9500])
Eg=7.78*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=1*10^-2; %% [cm^2]-> Area rivelatore

else if strcmp(S,'Diamond')                 

figure
plot(data(:,1)*10^9,smooth(data(:,2)*dB,10));
Eg=13*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=2.5*10^-2;
end
end
%% Insert the time shift and the other parameters
dt=8592; % [ns]

%%FLIGHT PATH%%
r=1.135; % [m]

c2=9*10^(16); % [m2/s2] c^2
c=3*10^8; % [m/s]
k=8.6*10^-5; % [eV] Cost Boltzmann
%%Type of particle
mp= 1*((938.27)/(9*10^16)); % [MeV/c2] col 4 sono alfa

e=1.6*10^-19; 


%%Flight path in cm
r_cm=r*100;

% Detector solid angle
Ang=S_d/(r_cm^2); 

Dr=0.01; %% uncertainty on the distance in cm
Dt=10^-(10); % uncertainty on time in s
dAng=2*S_d*(Dr)/(r_cm)^3

tl=r/c;% time of the light in sec

D_1=data(:,1)-(dt*10^-9)+tl; %% total time shift [sec]
D_2=(data(:,2)*dB) %% ampiezza segnale
 
% Select the "time" range of interest from TOF signal using the brush tool and call the variable "proton"
figure;
plot(D_1*10^9,smooth(D_2,10));
xlim([-200 400])
%% Load the parameters of the fits needed to convert the TOF signal in energy spectrum
%%Residual energy after the filter PROTONS:
param_filter=load('/Users/Ferdi/Ferdinanda Dati Prova/DatiProva Spettro/ConFiltroCheSiFermano/Res_6.5_Al_0.570_4MeV_protons.txt');

%% For particle which after the filter stop inside the detector active layer
%%select in the plot and call the variable proton
time_min_stop=proton(1,1); % ns
Ene_time_min_stop=(0.5*mp*r^2)/((time_min_stop*10^-9)^2) %1/2 mv2

Ene_time_max_stop=0.5;  % energia minima per cui le particelle passano il filtro e non si fermano (dentro il filtro)
time_max_stop=sqrt((0.5*mp*r^2)/(Ene_time_max_stop))*10^9  % [ns]


%time_max_stop=140;%ns
%Ene_time_max_stop=(0.5*mp*r^2)/((time_min_stop*10^-9)^2)

D0_2(:,1)=proton(:,1);
D0_2(:,2)=proton(:,2);

t_stop=D0_2(D0_2(:,1)>time_min_stop&D0_2(:,1)<time_max_stop,1)*10^-9;
V_stop=D0_2(D0_2(:,1)>time_min_stop&D0_2(:,1)<time_max_stop,2);
figure;
plot(t_stop,V_stop);
R=50;

Dr=1;
for(i=1:length(t_stop))
 beta_2(i)=(r/(t_stop(i)))/c;
 gamma_2(i)=1/(sqrt(1-beta_2(i).^2));
 E_inc_stop(i)=(gamma_2(i)-1)*mp*c^2;
%E_stop(i)=((mp*r*r)/(2*t_stop_2(i)*t_stop_2(i)));
sigma_inc2(i)=sqrt(((mp*r/t_stop(i))^2*(Dr^2))+((((mp*r^2)/(t_stop(i).^3)).^2)*(Dt)^2));
end;
figure
plot(E_inc_stop,V_stop)


%% Retrieve the residual energy after the filter using the parameters obtained from the fit->Poly9->param_filter
Eres_stop=param_filter(1)*(E_inc_stop.^9)+param_filter(2)*(E_inc_stop.^8)+param_filter(3)*(E_inc_stop.^7)+param_filter(4)*(E_inc_stop.^6)+param_filter(5)*(E_inc_stop.^5)+param_filter(6)*(E_inc_stop.^4)+param_filter(7)*(E_inc_stop.^3)+param_filter(8)*(E_inc_stop.^2)+param_filter(9)*(E_inc_stop)+param_filter(10);

%%
figure;
hold on
plot(E_inc_stop,Eres_stop)
%%
sigma_Eres=0.1;%Error in residual energy calculation-> straggling is about 100 KeV

for(i=1:length(t_stop))
dN_dE_stop(i)=((Eg*V_stop(i))/(e*R*Eres_stop(i)*Eres_stop(i)))*(0.5*t_stop(i)-Dt)
sigma_dN_stop(i)=(((Eg*V_stop(i))/(e*R*Ang))*sqrt((((t_stop(i))/(Eres_stop(i).^3)).^2*(sigma_Eres)^2)+((1/(Eres_stop(i).^2*2)).^2*Dt^2)+(t_stop(i)/(2*Eres_stop(i)*Ang)).^2*dAng^2));
err_per_stop(i)=(sigma_dN_stop(i)/dN_dE_stop(i))*100;
end;

%%
figure;
plot(Eres_stop,dN_dE_stop/Ang)
xlabel('Released proton energy [MeV] ','FontSize',20);
ylabel('dN/dE_{res} \Omega [MeV^{-1}sr^{-1}]','FontSize',20);

% Number
NN_stop=trapz(Eres_stop,dN_dE_stop)

%% Rebinning to have dN_dE inc vs Einc
sigmatot_2=0;
dEA_2=diff(Eres_stop);
for (i=1:length(Eres_stop)-1)
    dIntA_2(i)=dN_dE_stop(i)*dEA_2(i);
    sigma_2(i)=-(sigma_dN_stop(i)*dEA_2(i));
    sigmatot_2=sigma_2(i)+sigmatot_2;
end
%%Interval of incident energy
dEA_inc_2=diff(E_inc_stop);
%%dN/dE inc
for (i=1:length(E_inc_stop)-1)
    dN_dE_inc_stop(i)=dIntA_2(i)/dEA_inc_2(i);
    sigma_inc_2(i)=sigma_2(i)/dEA_inc_2(i);
end

Number_particles_omega_2=-trapz(E_inc_stop(1:length(E_inc_stop)-1),dN_dE_inc_stop)/Ang
Err_sist_rel_2=(sigmatot_2/Number_particles_omega_2)*100
Err_pois_2=sqrt(Number_particles_omega_2)
Err_tot_2=sqrt(sigmatot_2^2+Err_pois_2^2)
Err_rel_tot_2=(Err_tot_2/Number_particles_omega_2)*100
% Energy spectrum
figure
hold on;
plot(E_inc_stop(1:length(E_inc_stop)-1),dN_dE_inc_stop/Ang,'r')
xlabel('Incident proton energy [MeV] ','FontSize',20);
ylabel('dN/dE \Omega [MeV^{-1}sr^{-1}]','FontSize',20);
%%
clear
