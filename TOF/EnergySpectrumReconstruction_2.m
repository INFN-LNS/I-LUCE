%% Experimental data 
% C'è filtro e c'è lo strato morto e i protoni non si fermano dentro al rivelatore

% Rivelatore: SiC 
% Spessore Rivelatore:10 um
% Diametro Rivelatore= 
% Superficie Detector= 0.8 mm^2
% Distanza target-rivelatore= 850 cm

% Spessore filtro alluminio= 0

data = textread('C:\Users\Ferdi\OneDrive\Documenti\MATLAB\PALS2020\DATI\02112020\56004\C1scope2_5600400000.txt', '', 'delimiter', ',', ... 
              'headerlines',5,'emptyvalue', NaN);

% Commenta SiC o Diamante e controlla le aree S_d
           %S='SiC';
           S='Diamond'; 
           dB=1; %%-> Attenuation 20 DB->10; 10 DB->3; 6 DB->2; 13 DB-> 4.45
           
if  strcmp(S,'SiC')
   
figure;
plot(data(:,1)*10^9,smooth(data(:,2)*dB,5));
xlim([0 600])
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
Eg=7.78*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=0.8*10^-2; %% [cm^2]-> Display area rivelatore

else if strcmp(S,'Diamond')     
                    
figure;
plot(data(:,1)*10^9,smooth(data(:,2)*dB,1));
%xlim([-50 400])
Eg=13*10^-6 %% [MeV] electron-hole pair energy creation 13 eV for Diamond, 7.78 for SiC 
S_d=1.0*10^-2;
end
end

%% Insert the time shift and the other parameters
dt=195.5; % [ns]

%%FLIGHT PATH%%
r=0.85; % [m]

c2=9*10^(16); % [m2/s2] c^2
c=3*10^8; % [m/s]
k=8.6*10^-5; % [eV] Cost Boltzmann
%%Type of particle
mp= 1*((938.27)/(9*10^16)); % [MeV/c2] col 4 sono alfa

e=1.6*10^-19; 
%SiC
%Diamond %% 
%%Flight path in cm
r_cm=r*100;

%detector solid angle
Ang=S_d/(r_cm^2); 

Dr=0.01; %% uncertainty on the distance in cm
Dt=10^-(10); % uncertainty on time in s
dAng=2*S_d*(Dr)/(r_cm)^3

tl=r/c;% time of the light in sec

D_1=data(:,1)-(dt*10^-9)+tl; %%total time shift [sec]
D_2=(data(:,2)*dB)
 
% Select the "time" range of interest from TOF signal using the brush tool and call the variable "proton"
figure;
plot(D_1*10^9,smooth(D_2,5));
xlim([-100 400])
xlabel('TOF [ns] ','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
%% Load the parameters of the fits needed to convert the TOF signal in energy spectrum

param_eloss=load('\Users\Ferdi\OneDrive\Documenti\MATLAB\TESI\Dati per FIT\Senza filtro ma Strato Morto e NON si fermano\Res_10umSiC+SM_1.1_2.4MeV_protons.txt');


%% ENERGY SPECTRUM FOR PROTONS WHICH ALSO AFTER THE FILTER DON'T STOP WITHIN THE DETECTOR
time_min_loss=proton(1,1);% [ns] % Time corresponding to the maximum incident energy-> i.e. the first point of the variable proton
Ene_time_min_loss=(0.5*mp*r^2)/((time_min_loss*10^-9)^2)

%time_max_loss=68.5 % ns 
%Ene_time_max_loss=(0.5*mp*r^2)/((time_max_loss*10^-9)^2)
Ene_time_max_loss=1.1; % energia minima per cui le particelle passano lo strato morto e non si fermano
time_max_loss=sqrt((0.5*mp*r^2)/(Ene_time_max_loss))*10^9


%% Time corresponding to the "punch through" point-> particles with an energy less than this value 
% after the filter stop within the detector
D0_1(:,1)=proton(:,1);
D0_1(:,2)=proton(:,2);
%%select the interval of interest
t_loss=D0_1(D0_1(:,1)>time_min_loss&D0_1(:,1)<time_max_loss,1)*10^-9;
V_loss=D0_1(D0_1(:,1)>time_min_loss&D0_1(:,1)<time_max_loss,2);
figure;
plot(t_loss,V_loss);
Dt=10^-10; % time sampling
R=50; %resistance

%Incident energy in this energy range
for(i=1:length(t_loss))
 beta_1(i)=(r/(t_loss(i)))/c;
 gamma_1(i)=1/(sqrt(1-beta_1(i).^2));
 E_inc_loss(i)=(gamma_1(i)-1)*mp*c^2;
 sigma_inc1(i)=sqrt(((mp*r/t_loss(i))^2*(Dr^2))+((((mp*r^2)/(t_loss(i).^3)).^2)*(Dt)^2));
end;
figure;
plot(E_inc_loss,V_loss)
%% Retrieve the residual energy after the filter using the parameters obtained from the fit->Poly9->param_filter
%Eres_loss=param_filter(1)*(E_inc_loss.^9)+param_filter(2)*(E_inc_loss.^8)+param_filter(3)*(E_inc_loss.^7)+param_filter(4)*(E_inc_loss.^6)+param_filter(5)*(E_inc_loss.^5)+param_filter(6)*(E_inc_loss.^4)+param_filter(7)*(E_inc_loss.^3)+param_filter(8)*(E_inc_loss.^2)+param_filter(9)*(E_inc_loss)+param_filter(10);
Eres_loss=E_inc_loss; % come fosse caso senza filtro
%% To plot incident energy vs residual energy after the filter
figure;
plot(E_inc_loss,Eres_loss) % energia residua dopo il filtro
xlabel('Energy [MeV] ','FontSize',20);
ylabel('Residual energy [MeV]','FontSize',20);

%% Retrieve the energy loss within the detector using the parameter obtained from the fit-> param_eloss
%%if the fit is poly9
DE=param_eloss(1)*Eres_loss.^9+param_eloss(2)*Eres_loss.^8+param_eloss(3)*Eres_loss.^7+param_eloss(4)*Eres_loss.^6+param_eloss(5)*Eres_loss.^5+param_eloss(6)*Eres_loss.^4+param_eloss(7)*Eres_loss.^3+param_eloss(8)*Eres_loss.^2+param_eloss(9)*Eres_loss+param_eloss(10);

%% To plot incident energy vs energy loss in detector thickness
figure;
plot(E_inc_loss,DE,'.');
xlabel('Incident energy [MeV] ','FontSize',20);
ylabel('Energy loss in detector [MeV]','FontSize',20);
%% Reconstruct the energy spectrum
dDeltaS=0
sigma_E=0.1; %100 KeV-->straggling in the energy loss
%%Energy spectrum per solid angle
for(i=1:length(t_loss))
dN_dE_loss(i)=((Eg*V_loss(i))/(e*R*DE(i)*DE(i)))*(0.5*t_loss(i)-Dt)/Ang;
sigma_dN_loss(i)=(((Eg*V_loss(i))/(e*R*Ang))*sqrt((((t_loss(i))/(DE(i).^3)).^2*(sigma_E)^2)+((1/(DE(i).^2*2)).^2*Dt^2)+(t_loss(i)/(2*DE(i)*Ang)).^2*dAng^2));
err_per_loss(i)=(sigma_dN_loss(i)/dN_dE_loss(i))*100; % errore percentuale
end;

%%
%Number of protons in this energy range
NN_loss=trapz(DE,dN_dE_loss)
%clearfigure;
plot(DE,dN_dE_loss,'.');
%% Re-binning to have dN_dE_inc vs E_inc PLOT FINALE 
dEA_1=diff(DE);
sigmatot_1=0;
for (i=1:length(DE)-1)
    dIntA_1(i)=dN_dE_loss(i)*dEA_1(i);
    sigma_1(i)=-(sigma_dN_loss(i)*dEA_1(i));
    sigmatot_1=sigma_1(i)+sigmatot_1
    
end
%%Interval of incident energy
dEA_inc_1=-diff(E_inc_loss);
%%dN/dE inc
for (i=1:length(E_inc_loss)-1)
    dN_dE_inc_loss(i)=dIntA_1(i)/dEA_inc_1(i);
    sigma_inc_loss(i)=sigma_1(i)/dEA_inc_1(i);
end
Number_particles_omega_1=trapz(E_inc_loss(1:length(E_inc_loss)-1),dN_dE_inc_loss)
Err_sist_rel_1=(sigmatot_1/Number_particles_omega_1)*100
Err_pois_1=sqrt(abs(Number_particles_omega_1))
Err_tot_1=sqrt(sigmatot_1^2+Err_pois_1^2)
Err_rel_tot_1=(Err_tot_1/Number_particles_omega_1)*100

% Energy spectrum
figure
hold on
plot(E_inc_loss(1:length(E_inc_loss)-1),dN_dE_inc_loss,'m')
xlabel('Incident proton energy [MeV] ','FontSize',20);
ylabel('dN/dE \Omega [MeV^{-1}sr^{-1}]','FontSize',20);
%% With uncertainties
figure;
errorbar(E_inc_loss(1:length(E_inc_loss)-1),dN_dE_inc_1,sigma_inc_1)

%% For particle which after the filter stop inside the detector active layer
%%select in the plot and call the variable proton
time_min_stop=proton(1,1); % ns
Ene_time_min_stop=(0.5*mp*r^2)/((time_min_stop*10^-9)^2)

Ene_time_max_stop=2.8;
time_max_stop=sqrt((0.5*mp*r^2)/(Ene_time_max_stop))*10^9    


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
%Eres_stop=E_inc_stop;
%figure;
%plot(E_inc_stop,Eres_stop)

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
