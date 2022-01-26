%% Read signals from oscilloscope 1
 
TOF_Scope1=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/Scope1/Scope1_56028.csv', '', 'delimiter', ',', ... 
                'headerlines',0,'emptyvalue', NaN);
            
   Ch1_1(:,1)=TOF_Scope1(:,1);
   Ch1_1(:,2)=TOF_Scope1(:,2);
   
   Ch2_1(:,1)=TOF_Scope1(:,1);
   Ch2_1(:,2)=TOF_Scope1(:,3);
   
   Ch3_1(:,1)=TOF_Scope1(:,1);
   Ch3_1(:,2)=TOF_Scope1(:,4);
   
   Ch4_1(:,1)=TOF_Scope1(:,1);
   Ch4_1(:,2)=TOF_Scope1(:,5);
  
  %% This will plot the channels all together
% figure;
% plot(Ch1_1(:,1)*10^9,Ch1_1(:,2)*2,Ch2_1(:,1)*10^9,Ch2_1(:,2),Ch3_1(:,1)*10^9,Ch3_1(:,2),Ch4_1(:,1)*10^9,Ch4_1(:,2));
% xlabel('TOF [ns]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)
%% Ch1 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
% figure;
% plot(Ch1_1(:,1)*10^9,smooth(Ch1_1(:,2),30))
%%
% dt_1_1=235.7;%ns   <---- time shift of the photopeak
% d1_1=0.85; %%%det to Target distance in m
% Att1_1=10 ; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
% tl_1_1=(d1_1/(3*10^8))*10^9 %light travel time
% Ch1_1(:,2)=Ch1_1(:,2)*Att1_1;
% time_Ch1_1=Ch1_1(:,1)*10^9-dt_1_1+tl_1_1 %ns TOF corrected for the time of light
% 
% % energy plot -> 
% %select with brush tool a tof spectrum region (do not select the photopeak)
% %to convert the tof in energy, save the variable as "proton_1_1"
% figure;
% plot(time_Ch1_1,smooth(Ch1_1(:,2),30),'m')
% xlabel('TOF [ns]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
% mP=((938.27)/(9*10^16));
% c=3*10^8;
% S1_1=0.8*10^-6; %detector's surface
%  for(i=1:length(proton_1_1))
% betaP_1_1(i)=(d1_1/(proton_1_1(i,1)*10^-9))/c;
% gammaP_1_1(i)=1/(sqrt(1-betaP_1_1(i).^2));
% %%RELATIVISTICA
% EP_1_1(i)=(gammaP_1_1(i)-1)*mP*c^2;
% I_1_1(i)=proton_1_1(i,2)*(d1_1^2)/(50*S1_1);
%  end
%  
% figure;
% hold on
% plot(EP_1_1,smooth(proton_1_1(:,2),1));
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch2 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch2_1(:,1)*10^9,smooth(Ch2_1(:,2),5))
%%
dt_2_1=169.7; %ns   <---- time shift of the photopeak
d2_1=1.65;  %%%det to Target distance in m
tl_2_1=(d2_1/(3*10^8))*10^9 
Att2_1=10; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch2_1(:,2)=Ch2_1(:,2)*Att2_1;
time_Ch2_1=Ch2_1(:,1)*10^9-dt_2_1+tl_2_1
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_2_1"
figure;
plot(time_Ch2_1,smooth(Ch2_1(:,2),40))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S2_1=0.8*10^-6; %detector's surface
 for(i=1:length(proton_2_1))
betaP_2_1(i)=(d2_1/(proton_2_1(i,1)*10^-9))/c;
gammaP_2_1(i)=1/(sqrt(1-betaP_2_1(i).^2));
%%RELATIVISTICA
EP_2_1(i)=(gammaP_2_1(i)-1)*mP*c^2;
I_2_1(i)=proton_2_1(i,2)*(d2_1^2)/(50*S2_1);
 end
 
% figure;
% plot(EP_2_1,smooth(proton_2_1(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch3 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch3_1(:,1)*10^9,-smooth(Ch3_1(:,2),10))
%%
dt_3_1=311.5; %ns   <---- time shift of the photopeak
d3_1=1.62; %%%det to Target distance in m
tl_3_1=(d3_1/(3*10^8))*10^9
Att3_1=10; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch3_1(:,2)=-Ch3_1(:,2)*Att3_1;
time_Ch3_1=Ch3_1(:,1)*10^9-dt_3_1+tl_3_1
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_3_1"
figure;
%hold on;
plot(time_Ch3_1,smooth(Ch3_1(:,2),1))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S3_1=225*10^-6;
 for(i=1:length(proton_3_1))
betaP_3_1(i)=(d3_1/(proton_3_1(i,1)*10^-9))/c;
gammaP_3_1(i)=1/(sqrt(1-betaP_3_1(i).^2));
%%RELATIVISTICA
EP_3_1(i)=(gammaP_3_1(i)-1)*mP*c^2;
I_3_1(i)=proton_3_1(i,2)*(d3_1^2)/(50*S3_1);
 end
 
% figure;
% hold on
% plot(EP_3_1,smooth(proton_3_1(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch4-> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch4_1(:,1)*10^9,smooth(Ch4_1(:,2),10))
%%
dt_4_1=144.5; %ns   <---- time shift of the photopeak
d4_1=1.51 %%%det to Target distance in m
tl_4_1=(d4_1/(3*10^8))*10^9
Att4_1=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch4_1(:,2)=Ch4_1(:,2)*Att4_1;
time_Ch4_1=Ch4_1(:,1)*10^9-dt_4_1+tl_4_1
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_4_1"
figure;
plot(time_Ch4_1,smooth(Ch4_1(:,2),20))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S4_1=2.5*10^-6;
 for(i=1:length(proton_4_1))
betaP_4_1(i)=(d4_1/(proton_4_1(i,1)*10^-9))/c;
gammaP_4_1(i)=1/(sqrt(1-betaP_4_1(i).^2));
%%RELATIVISTICA
EP_4_1(i)=(gammaP_4_1(i)-1)*mP*c^2;
I_4_1(i)=proton_4_1(i,2)*(d4_1^2)/(50*S4_1);
 end
 
% figure;
% hold on
% plot(EP_4_1,smooth(proton_4_1(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

 %% Read signals from oscilloscope 2
           
Ch1_2=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/Scope2/C1Scope2_5602800000.txt', '', 'delimiter', ',', ... 
                'headerlines',5,'emptyvalue', NaN);
Ch2_2=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/Scope2/C2scope2_5602800000.txt', '', 'delimiter', ',', ... 
                'headerlines',5,'emptyvalue', NaN);
Ch3_2=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/Scope2/C3scope2_5602800000.txt', '', 'delimiter', ',', ... 
              'headerlines',5,'emptyvalue', NaN);
Ch4_2=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/Scope2/C4scope2_5602800000.txt', '', 'delimiter', ',', ... 
                'headerlines',5,'emptyvalue', NaN);
           
  
  %% This will plot the channels all together
% figure;
% plot(Ch1_2(:,1)*10^9,Ch1_2(:,2)*2,Ch2_2(:,1)*10^9,Ch2_2(:,2),Ch3_2(:,1)*10^9,Ch3_2(:,2),Ch4_2(:,1)*10^9,Ch4_2(:,2));
% xlabel('TOF [ns]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)
%% Ch1 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch1_2(:,1)*10^9,smooth(Ch1_2(:,2),10))
%%
dt_1_2=147.8;%ns   <---- time shift of the photopeak
d1_2=0.85; %%%det to Target distance in m
Att1_2=1 ; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
tl_1_2=(d1_2/(3*10^8))*10^9 %light travel time
Ch1_2(:,2)=Ch1_2(:,2)*Att1_2;
time_Ch1_2=Ch1_2(:,1)*10^9-dt_1_2+tl_1_2 %ns TOF corrected for the time of light

% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_1_2"
figure;
plot(time_Ch1_2,smooth(Ch1_2(:,2),10),'m')
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
 c=3*10^8;
S1_2=1*10^-6; %detector's surface
 for(i=1:length(proton_1_2))
betaP_1_2(i)=(d1_2/(proton_1_2(i,1)*10^-9))/c;
gammaP_1_2(i)=1/(sqrt(1-betaP_1_2(i).^2));
%%RELATIVISTICA
EP_1_2(i)=(gammaP_1_2(i)-1)*mP*c^2;
I_1_2(i)=proton_1_2(i,2)*(d1_2^2)/(50*S1_2);
 end
 
% figure;
% hold on
% plot(EP_1_2,smooth(proton_1_2(:,2),1));
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch3 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch3_2(:,1)*10^9,smooth(Ch3_2(:,2),10))
%%
dt_3_2=121.6; %ns   <---- time shift of the photopeak
d3_2=1.76; %%%det to Target distance in m
tl_3_2=(d3_2/(3*10^8))*10^9
Att3_2=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch3_2(:,2)=Ch3_2(:,2)*Att3_2;
time_Ch3_2=Ch3_2(:,1)*10^9-dt_3_2+tl_3_2
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_3_2"
figure;
%hold on;
plot(time_Ch3_2,smooth(Ch3_2(:,2),20))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S3_2=1*10^-6;
 for(i=1:length(proton_3_2))
betaP_3_2(i)=(d3_2/(proton_3_2(i,1)*10^-9))/c;
gammaP_3_2(i)=1/(sqrt(1-betaP_3_2(i).^2));
%%RELATIVISTICA
EP_3_2(i)=(gammaP_3_2(i)-1)*mP*c^2;
I_3_2(i)=proton_3_2(i,2)*(d3_2^2)/(50*S3_2);
 end
 
% figure;
% hold on
% plot(EP_3_2,smooth(proton_3_2(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch4-> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch4_2(:,1)*10^9,smooth(Ch4_2(:,2),10))
%%
dt_4_2=106.6; %ns   <---- time shift of the photopeak
d4_2=1.46 %%%det to Target distance in m
tl_4_2=(d4_2/(3*10^8))*10^9
Att4_2=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch4_2(:,2)=Ch4_2(:,2)*Att4_2;
time_Ch4_2=Ch4_2(:,1)*10^9-dt_4_2+tl_4_2
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_4_2"
figure;
plot(time_Ch4_2,smooth(Ch4_2(:,2),10))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S4_2=2.5*10^-6;
 for(i=1:length(proton_4_2))
betaP_4_2(i)=(d4_2/(proton_4_2(i,1)*10^-9))/c;
gammaP_4_2(i)=1/(sqrt(1-betaP_4_2(i).^2));
%%RELATIVISTICA
EP_4_2(i)=(gammaP_4_2(i)-1)*mP*c^2;
I_4_2(i)=proton_4_2(i,2)*(d4_2^2)/(50*S4_2);
 end
 
% figure;
% hold on
% plot(EP_4_2,smooth(proton_4_2(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

 %%  Read signals from oscilloscope 3
 
 TOF_Scope3=textread('/Users/Beatrice/Desktop/UNIVERSITA/MAGISTRALE/TESI MAGISTRALE/analisi dati/TOF/TOF_DATA/05112020/scope3/scope3_56028.csv', '', 'delimiter', ',', ... 
                'headerlines',22,'emptyvalue', NaN);
            
   Ch1_3(:,1)=TOF_Scope3(:,1);
   Ch1_3(:,2)=TOF_Scope3(:,2);
   
   Ch2_3(:,1)=TOF_Scope3(:,1);
   Ch2_3(:,2)=TOF_Scope3(:,3);
   
   Ch3_3(:,1)=TOF_Scope3(:,1);
   Ch3_3(:,2)=TOF_Scope3(:,4);
   
   Ch4_3(:,1)=TOF_Scope3(:,1);
   Ch4_3(:,2)=TOF_Scope3(:,5);
  
  %% This will plot the channels all together
% figure;
% plot(Ch1_3(:,1)*10^9,Ch1_3(:,2)*2,Ch2_3(:,1)*10^9,Ch2_3(:,2),Ch3_3(:,1)*10^9,Ch3_3(:,2),Ch4_3(:,1)*10^9,Ch4_3(:,2));
% xlabel('TOF [ns]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)
%% Ch1 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch1_3(:,1)*10^9,-smooth(Ch1_3(:,2),30))
%%
dt_1_3=102.8;%ns   <---- time shift of the photopeak
d1_3=1.495; %%%det to Target distance in m
Att1_3=1 ; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
tl_1_3=(d1_3/(3*10^8))*10^9 %light travel time
Ch1_3(:,2)=Ch1_3(:,2)*Att1_3;
time_Ch1_3=Ch1_3(:,1)*10^9-dt_1_3+tl_1_3 %ns TOF corrected for the time of light

% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_1_3"
figure;
plot(time_Ch1_3,-smooth(Ch1_3(:,2),10),'m')
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S1_3=4*10^-6;
 for(i=1:length(proton_1_3))
betaP_1_3(i)=(d1_3/(proton_1_3(i,1)*10^-9))/c;
gammaP_1_3(i)=1/(sqrt(1-betaP_1_3(i).^2));
%%RELATIVISTICA
EP_1_3(i)=(gammaP_1_3(i)-1)*mP*c^2;
I_1_3(i)=proton_1_3(i,2)*(d1_3^2)/(50*S1_3);
 end
 
% figure;
% plot(EP_1_3,smooth(proton_1_3(:,2),1));
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch2 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch2_3(:,1)*10^9,smooth(Ch2_3(:,2),20))
%%
dt_2_3=130.6; %ns   <---- time shift of the photopeak
d2_3=1.65;  %%%det to Target distance in m
tl_2_3=(d2_3/(3*10^8))*10^9 
Att2_3=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch2_3(:,2)=Ch2_3(:,2)*Att2_3;
time_Ch2_3=Ch2_3(:,1)*10^9-dt_2_3+tl_2_3
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_2_3"
figure;
plot(time_Ch2_3,smooth(Ch2_3(:,2),20))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S2_3=3.14*10^-6;
 for(i=1:length(proton_2_3))
betaP_2_3(i)=(d2_3/(proton_2_3(i,1)*10^-9))/c;
gammaP_2_3(i)=1/(sqrt(1-betaP_2_3(i).^2));
%%RELATIVISTICA
EP_2_3(i)=(gammaP_2_3(i)-1)*mP*c^2;
I_2_3(i)=proton_2_3(i,2)*(d2_3^2)/(50*S2_3);
 end
 
% figure;
% plot(EP_2_3,smooth(proton_2_3(:,2),1)*1);
% xlabel('Proton energy [MeV]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch3 -> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
figure;
plot(Ch3_3(:,1)*10^9,smooth(Ch3_3(:,2),30))
%%
dt_3_3=124.6; %ns   <---- time shift of the photopeak
d3_3=0.80; %%%det to Target distance in m
tl_3_3=(d3_3/(3*10^8))*10^9
Att3_3=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
Ch3_3(:,2)=Ch3_3(:,2)*Att3_3;
time_Ch3_3=Ch3_3(:,1)*10^9-dt_3_3+tl_3_3
% energy plot -> 
%select with brush tool a tof spectrum region (do not select the photopeak)
%to convert the tof in energy, save the variable as "proton_3_3"
figure;
%hold on;
plot(time_Ch3_3,smooth(Ch3_3(:,2),30))
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
mP=((938.27)/(9*10^16));
c=3*10^8;
S3_3=16*10^-6;
for(i=1:length(proton_3_3))
betaP_3_3(i)=(d3_3/(proton_3_3(i,1)*10^-9))/c;
gammaP_3_3(i)=1/(sqrt(1-betaP_3_3(i).^2));
%%RELATIVISTICA
EP_3_3(i)=(gammaP_3_3(i)-1)*mP*c^2;
I_3_3(i)=proton_3_3(i,2)*(d3_3^2)/(50*S3_3);
end
 
% figure;
% %plot(EP_3_3,smooth(proton_3_3(:,2),1)*1);
% plot(EP_3_3,I_3_3);
% xlabel('Proton energy [MeV]','FontSize',20);
% %ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)

%% Ch4-> identify the time shift dt (ns) to shift photopeak to zero and put it as input in the next section
% figure;
% plot(Ch4_3(:,1)*10^9,-smooth(Ch4_3(:,2),10))
%%
% dt_4_3=122.4; %ns   <---- time shift of the photopeak
% d4_3=1.495 %%%det to Target distance in m
% tl_4_3=(d4_3/(3*10^8))*10^9
% Att4_3=1; %% attenuation use Att=10->20 dB Att=2 ->6 dB Att=1.41-> 3 dB  Att=3.16-> dB=10
% Ch4_3(:,2)=-Ch4_3(:,2)*Att4_3;
% time_Ch4_3=Ch4_3(:,1)*10^9-dt_4_3+tl_4_3
% % energy plot -> 
% %select with brush tool a tof spectrum region (do not select the photopeak)
% %to convert the tof in energy, save the variable as "proton_4_3"
% figure;
% plot(time_Ch4_3,smooth(Ch4_3(:,2),10))
% xlabel('TOF [ns]','FontSize',20);
% ylabel('Amplitude [Volt]','FontSize',20);
% set(gca,'fontsize',20)
%% This convert the TOF values in proton kinetic energy 
% mP=((938.27)/(9*10^16));
% c=3*10^8;
% %S4_3=1*10^-6;
% for(i=1:length(proton_4_3))
% betaP_4_3(i)=(d4_3/(proton_4_3(i,1)*10^-9))/c;
% gammaP_4_3(i)=1/(sqrt(1-betaP_4_3(i).^2));
% %%RELATIVISTICA
% EP_4_3(i)=(gammaP_4_3(i)-1)*mP*c^2;
% %I_4_3(i)=proton_4_3(i,2)*(d4_3^2)/(50*S4_3);
% end
%  
% % figure;
% % plot(EP_4_3,smooth(proton_4_3(:,2),1)*1);
% % xlabel('Proton energy [MeV]','FontSize',20);
% % ylabel('Amplitude [Volt]','FontSize',20);
% % set(gca,'fontsize',20)

%% PLOT ALL THE TOF SIGNALS TOGETHER
figure;
plot(time_Ch1_1,Ch1_1(:,2),'DisplayName', 'S1 0° b')
hold on
plot(time_Ch2_1,Ch2_1(:,2),'DisplayName', 'S2 30° f')
%plot(time_Ch3_1,Ch3_1(:,2),'DisplayName', '')
%plot(time_Ch4_1,smooth(Ch4_1(:,2),20),'DisplayName', 'D5 60° f')

plot(time_Ch1_2,Ch1_2(:,2),'DisplayName', 'D7 30° b')
%plot(time_Ch3_2,Ch3_2(:,2),'DisplayName', 'D3 60° f')
plot(time_Ch4_2,smooth(Ch4_2(:,2),20),'DisplayName', 'D4 45° f')

plot(time_Ch1_3,Ch1_3(:,2),'DisplayName', 'D1 0° f')
plot(time_Ch2_3,Ch2_3(:,2),'DisplayName', 'D2 30° f')
plot(time_Ch3_3,Ch3_3(:,2),'DisplayName', 'D6 30° b')
%plot(time_Ch4_3,smooth(Ch4_3(:,2),20),'DisplayName', 'Si0 0° f')
hold off
xlabel('TOF [ns]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)
%%  PLOT ALL THE TOF SIGNALS  TOGETHER NORMALIZING THE TOF AT 1M 
figure;
plot(time_Ch1_1/d1_1,Ch1_1(:,2),'DisplayName', 'S1 0° b')
hold on
plot(time_Ch2_1/d2_1,Ch2_1(:,2),'DisplayName', 'S2 30° f')
%plot(time_Ch3_1/d3_1,Ch3_1(:,2),'DisplayName','')
%plot(time_Ch4_1/d4_1,smooth(Ch4_1(:,2),20),'DisplayName', 'D5 60° f')

plot(time_Ch1_2/d1_2,Ch1_2(:,2),'DisplayName', 'D7 30° b')
%plot(time_Ch3_2/d3_2,Ch3_2(:,2),'DisplayName', 'D3 60° f')
plot(time_Ch4_2/d4_2,smooth(Ch4_2(:,2),20),'DisplayName', 'D4 45° f')

plot(time_Ch1_3/d1_3,Ch1_3(:,2),'DisplayName', 'D1 0° f')
plot(time_Ch2_3/d2_3,Ch2_3(:,2),'DisplayName', 'D2 30° f')
plot(time_Ch3_3/d3_3,Ch3_3(:,2),'DisplayName', 'D6 30° b')
%plot(time_Ch4_3/d4_3,smooth(Ch4_3(:,2),20),'DisplayName', 'Si0 0° f')
hold off

%% PLOT ALL THE PROTON ENERGIES TOGETHER
figure;
%plot(EP_1_1,smooth(proton_1_1(:,2),1),'DisplayName', 'S1 0° b')
plot(EP_2_1,smooth(proton_2_1(:,2),1),'DisplayName', 'S2 30° f')
hold on
plot(EP_3_1,smooth(proton_3_1(:,2),1),'DisplayName','D_ poli 0° b')
plot(EP_4_1,smooth(proton_4_1(:,2),1),'DisplayName', 'D5 60° f')

plot(EP_1_2,smooth(proton_1_2(:,2),1),'DisplayName', 'D7 30° b')
plot(EP_3_2,smooth(proton_3_2(:,2),1),'DisplayName', 'D3 60° f')
plot(EP_4_2,smooth(proton_4_2(:,2),1),'DisplayName', 'D4 45° f')

plot(EP_1_3,smooth(proton_1_3(:,2),1),'DisplayName', 'D1 0° f')
plot(EP_2_3,smooth(proton_2_3(:,2),1),'DisplayName', 'D2 30° f')
plot(EP_3_3,smooth(proton_3_3(:,2),1),'DisplayName', 'D6 30° b')
%plot(EP_4_3,smooth(proton_4_3(:,2),1),'DisplayName', 'Si0 0° f')
hold off
xlabel('Proton energy [MeV]','FontSize',20);
ylabel('Amplitude [Volt]','FontSize',20);
set(gca,'fontsize',20)

%% PLOT ALL THE NORMALIZED PROTON ENERGIES TOGETHER
figure;
%plot(EP_1_1,I_1_1,'DisplayName', 'S1 0° b')
plot(EP_2_1,I_2_1,'DisplayName', 'S2 30° f')
hold on
plot(EP_3_1,I_3_1,'DisplayName','D_ poli 0° b')
plot(EP_4_1,I_4_1,'DisplayName', 'D5 60° f')

plot(EP_1_2,I_1_2,'DisplayName', 'D7 30° b')
plot(EP_3_2,I_3_2,'DisplayName', 'D3 60° f')
plot(EP_4_2,I_4_2,'DisplayName', 'D4 45° f')

plot(EP_1_3,I_1_3,'DisplayName', 'D1 0° f')
plot(EP_2_3,I_2_3,'DisplayName', 'D2 30° f')
plot(EP_3_3,I_3_3,'DisplayName', 'D6 30° b')
%plot(EP_4_3,I_4_3,'DisplayName', 'Si0 0° f')
hold off
xlabel('Proton energy [MeV]','FontSize',20);
ylabel('Normalized Current [A/sr]','FontSize',20);
set(gca,'fontsize',20)
 %% Delete all the variables at the end of the analysis 
clear proton_1_1;
clear EP_1_1;
clear I_1_1;
clear gammaP_1_1; 
clear betaP_1_1;

clear proton_2_1;
clear EP_2_1;
clear I_2_1;
clear gammaP_2_1;
clear betaP_2_1;

clear proton_3_1;
clear EP_3_1;
clear I_3_1;
clear gammaP_3_1;
clear betaP_3_1;

clear proton_4_1;
clear EP_4_1;
clear I_4_1;
clear gammaP_4_1;
clear betaP_4_1;

clear proton_1_2;
clear EP_1_2;
clear I_1_2;
clear gammaP_1_2;
clear betaP_1_2;

clear proton_3_2;
clear EP_3_2;
clear I_3_2;
clear gammaP_3_2;
clear betaP_3_2;

clear proton_4_2;
clear EP_4_2;
clear I_4_2;
clear gammaP_4_2;
clear betaP_4_2;

clear proton_1_3;
clear EP_1_3;
clear I_1_3;
clear gammaP_1_3;
clear betaP_1_3;

clear proton_2_3;
clear EP_2_3;
clear I_2_3;
clear gammaP_2_3; 
clear betaP_2_3;

clear proton_3_3;
clear EP_3_3;
clear I_3_3;
clear gammaP_3_3; 
clear betaP_3_3;

clear proton_4_3;
clear EP_4_3;
clear I_4_3;
clear gammaP_4_3;
clear betaP_4_3;