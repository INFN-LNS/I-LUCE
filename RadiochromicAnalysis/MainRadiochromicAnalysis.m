clear all;
close all;

%% Carico il file con la matrice degli stopping power

StoppingPowerMatrix=load('StoppingPowerMatrix_64MeV.mat'); 

%% Carico la distribuzione di dose sperimentale
%

% DatiMarkus=load('Markus_FE_21_02_2019.txt');  %60MeV
DatiMarkus = load('Markus_FE_30MeV1.txt'); %30MeV

%%
% Pongo uguale a zero tutti i valori negativi della terza colonna
% del file della distribuzione di dose. La terza colonna di questo file 
% contiene i valori normalizzati della dose 
%
[NegativeIndex] = find(DatiMarkus(:,3) < 0 );
DatiMarkus(NegativeIndex,3) = 0; 

%%
% Carico una spettro di energia dei protoni generato da una simulazione 
% Monte Carlo.  
% This will be used to compare the calculated final spectra
%
load('ProtonSpectraEntrance_60MeV.mat');

%% Conversion factors and common parameters
%
ScannerNominalResolution = 150; % in dpi (pixels per inch)
PixelDimension = 2.54/ScannerNominalResolution; % one side of the pixel in cm
PixelSurface = PixelDimension^2;

% Calculation of the beam spot surface in cm2
%
BeamDiameter = 1.5; % in cm (if circular....)
BeamSpotSurface = pi* (BeamDiameter/2)^2;

% Beam spot surface in number of Pixels
%
BeamSpotSurfaceInPixels = BeamSpotSurface/PixelSurface;

%Spot_px=BeamSpotSurface/(PxcmFactor150^2); % sezione trasversale del fascio in #pixel usando un collimatore con diametro da 15 mm
sensibility_intensity=1; % ? la sensibilit? di misura con cui valuto i livelli di intensit?

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Carico le immagini che servono per la calibrazione, il fondo e lo
%%%%%%% stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('Seleziona una immagine qualsiasi nella cartella della calibrazione')
% [FileName_Calibration,PathName_Calibration]=uigetfile('*.tif','Seleziona una immagine qualsiasi nella cartella della calibrazione');
% TotCal=input('\n Inserisci il numero di file che vuoi usare per effettuare la calibrazione\n'); %%numero totale di immagini della calibrazione

rect_calibration=[77.5100 87.5100 23.9800 22.9800]; % FE 60MeV
Numero_di_pixel_ROI_Cal=round(rect_calibration(1,3))*(rect_calibration(1,4)); %%round arrotonda a interi ==> la roi ? quadrata
AreaRectCal_cm2=Numero_di_pixel_ROI_Cal*PixelSurface; %area della ROI in cm2 
% 
% 
%%
disp('Seleziona immagine da usare come fondo')
[FileName_Fondo,PathName_Fondo]=uigetfile('*.tif','Seleziona immagine da usare come fondo');
fullchosenfile_Fondo = [PathName_Fondo FileName_Fondo];
Fondo=imread(fullchosenfile_Fondo);


disp('Seleziona una immagine qualsiasi nella cartella dei dati che vuoi analizzare')
[FileName_Stack,PathName_Stack]=uigetfile('*.tif','Seleziona una immagine qualsiasi nella cartella dei dati che vuoi analizzare.');
TotStack=input('\n Inserisci il numero di gaf dello stack\n');
rect_stack=[78.5100   88.5100   22.9800   20.9800]; %FE 60 MeV  
Numero_di_pixel_ROI_Stack=round(rect_stack(1,3))*round(rect_stack(1,4));
AreaRectStack=[round(rect_stack(1,3))*PixelDimension*round(rect_stack(1,4))*PixelDimension]; %area della ROI in cm2 --> MODO 1

base=round(rect_stack(1,3))*PixelDimension; % cm
altezza=round(rect_stack(1,4))*PixelDimension; % cm
Err_AreaRectStack=sqrt((base)^2+(altezza)^2)*PixelDimension; % errore della misura dell'area della ROI in cm2

%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%% Calcolo l'intensita' media in ogni ROI di ogni gaf della
% % %%%%%%% calibrazione
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % for i=1:TotCal
% %        if i<10 %se le ultime due cifre del nome del file sono <10  
% %           AT=int2str(i); %trasforma la i in una stringa
% %           A=strcat('0',AT); %Qui metto lo zero davanti
% %        else 
% %           A=int2str(i);
% %        end
% % 
% % FileName_Calibration=strcat('EBT30',A,'.tif'); %S=EBT30 - A=i - Type=.tif
% % fullchosenfile = [PathName_Calibration FileName_Calibration];
% % GAFCAL=imread(fullchosenfile);
% % 
% % 
Xmin=round(rect_calibration(1,1));
Xmax=round(rect_calibration(1,1)+rect_calibration(1,3));
Ymin=round(rect_calibration(1,2));
Ymax=round(rect_calibration(1,2)+rect_calibration(1,4));
% % 
% % %j,k servono per scorrere lungo la ROI (secondo una serpentina) t serve per passare al pixel successivo
    t=1;
    for j=Xmin:Xmax
        for k=Ymin:Ymax
         
            %Vettore_Intensita_Gaf_Calibrazione(t,1)=GAFCAL(round(j),round(k),1); %pixelvaluered=PVR %1 ? la matrice che conserva il rosso
            Vettore_Intensita_Fondo(t,1)=Fondo(round(j),round(k),1);
            t=t+1;
            
        end
    end
% %                
% %    Average_Gaf_Cal(i,1)=mean(Vettore_Intensita_Gaf_Calibrazione); 
% %    Dev_Standard_Gaf_Cal(i,1)=std(Average_Gaf_Cal)./sqrt(Numero_di_pixel_ROI_Cal);
% %     
% % end
%%
Average_Fondo=mean(Vettore_Intensita_Fondo);
Double_Vettore_Intensita_Fondo=double(Vettore_Intensita_Fondo);
Dev_standard_Fondo=std(Double_Vettore_Intensita_Fondo)./sqrt(Numero_di_pixel_ROI_Cal);


%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%% Converto i valori in densita' ottica e calcolo l'errore
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % DensitaOttica_Calibrazione=log10(Average_Fondo./Average_Gaf_Cal);
% % Dose=[55.45 105.08 205.73 404.95 805.73 1004.86 1505.07 0];   % in cGy - FE - 60MeV
% % Err_Dose=3*Dose./100; 
% % 
% % syms x y
% % f=log10(x/y);
% % dfx=diff(f,x); % derivata parziale rispetto a x
% % dfy=diff(f,y); % derivata parziale rispetto a y
% % 
% % Derivata_parziale_x=vpa(subs(dfx,x,Average_Fondo)); % calcolo il valore numerico per l'unico valore di x
% % Derivata_parziale_y=vpa(subs(dfy,y,Average_Gaf_Cal));
% % 
% % Err_densitaottica_calibrazione=sqrt((Derivata_parziale_x.^2).*((Dev_standard_Fondo+sensibility_intensity).^2)+(Derivata_parziale_y.^2).*((Dev_Standard_Gaf_Cal+sensibility_intensity).^2)); 
% % Err_lettura_scanner=2*DensitaOttica_Calibrazione./100; % calcolo l'errore del 2% sulla lettura dello scanner

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CApplico la curva di calibrazione
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
a=575.7;
b=2.396; 
c=-1557;
d=0.4782;
e=981.7;

Funzione_FIT_DoppioEsponenziale=@(x) a*exp(x*b)+c*exp(x*d)+e; %la @x si usa per dirgli quale e' la variabile % QUI APPLICHIAMO LA CALIBRAZiONE

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calcolo l'intensit? media in ogni ROI di ogni gaf della
%%%%%%% stack e costruisco la matrice StackParameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StackParameter=zeros(TotStack,8);
GafParameter=zeros(TotStack,7);
DO_GAF=zeros(Numero_di_pixel_ROI_Stack,TotStack); % ? una matrice in cui ogni colonna ? relativa all'i-esimo GAF dello stack, 
                                               % mentre ogni riga contiene la DO corrispondente al singolo pixel della 
                                               % ROI usata per analizzare lo stack
Dev_Standard_DO_GAF=zeros(TotStack,1);

for i=1:TotStack 
    if i<10    
        
       A=strcat('0',int2str(i));
    else 
       A=int2str(i);
    end

   FileName_Stack=strcat('EBT30',A,'.tif');
   fullchosenfileStack = [PathName_Stack FileName_Stack];
   GAFSTACK=imread(fullchosenfileStack);
    

    Xmin=round(rect_stack(1,1));
    Xmax=round(rect_stack(1,1)+rect_stack(1,3));
    Ymin=round(rect_stack(1,2));
    Ymax=round(rect_stack(1,2)+rect_stack(1,4));
    
t=1;
m=1;

    for j=Xmin:Xmax
        for k=Ymin:Ymax
            PVR=double(GAFSTACK(j,k,1)); %pixelvaluered=PVR %1 ? la matrice che conserva il rosso
            
            
             Livello_Minimo=Average_Fondo;
            
            if PVR<Livello_Minimo %check per capire se ci sono macchie 
                
                %stabilisco il valore di fondo come livello minimo se il
                %pixel ha un valore minore del minimo (ovvero del fondo)
                %allora il valore minimo lo aggiorno
                
               Livello_Minimo=PVR;
              
               
               Matrice_Intensita_Gaf_Stack(t,i)=PVR;
               StackParameter(i,1)=StackParameter(i,1)+PVR; 
               StackParameter(i,2)=StackParameter(i,2)+1;
               t=t+1;
               
               
DO_PVR=log10(Average_Fondo/PVR); %qui calcoliamo la densita ottica del singolo pixel mentre in StackParameter(i,7) calcoliamo il valor medio
            DO_GAF(m,i)=DO_PVR;    
            DoseSinglePixel= Funzione_FIT_DoppioEsponenziale(DO_PVR); %%PR ? la curva di calibrazione quindi trovo la dose corrispondente alla densit? ottica di ogni singolo pixel
           
            
            GafParameter(i,1)=GafParameter(i,1)+DoseSinglePixel; %qui abbiamo la somma della dose di ogni pixel della roi
            GafParameter(i,2)=GafParameter(i,2)+1; %qui abbiamo il numero di pixel
m=m+1;


            else
                
               Matrice_Intensita_Gaf_Stack(t,i)=PVR; 
               StackParameter(i,1)=StackParameter(i,1)+PVR; %PVR per ogni ROI di ogni GAF
               StackParameter(i,2)=StackParameter(i,2)+1; %numero di pixel della ROI
t=t+1;
               

DO_PVR=log10(Average_Fondo/PVR); %qui calcoliamo la densita ottica del singolo pixel mentre in StackParameter(i,7) calcoliamo il valor medio
            DO_GAF(m,i)=DO_PVR;    
            DoseSinglePixel= Funzione_FIT_DoppioEsponenziale(DO_PVR); %%PR ? la curva di calibrazione quindi trovo la dose corrispondente alla densit? ottica di ogni singolo pixel
           
            
            GafParameter(i,1)=GafParameter(i,1)+DoseSinglePixel; %qui abbiamo la somma della dose di ogni pixel della roi
            GafParameter(i,2)=GafParameter(i,2)+1; %qui abbiamo il numero di pixel
m=m+1;

            end
        end
    end
    
    StackParameter(i,3)=Livello_Minimo;  %PVR del fondo
    StackParameter(i,4)=i; %numero del GAF
    StackParameter(i,5)=StackParameter(i,1)./StackParameter(i,2); %  valore medio dei PVR della ROI
    StackParameter(i,6)=std(Matrice_Intensita_Gaf_Stack(:,i))./sqrt(Numero_di_pixel_ROI_Stack); % deviazione standard sui PVR
    StackParameter(i,7)=log10(Average_Fondo/StackParameter(i,5)); %  DO media, ovvero calcolata sulla media dell'intensita della ROI 
    StackParameter(i,8)=Funzione_FIT_DoppioEsponenziale(StackParameter(i,7)); % la DOSE media in cGy calcolata considerando la DO media 
    
    
    GafParameter(i,3)=(GafParameter(i,1)*1*0.00001)*((PixelDimension)^2);  %% ATTENZIONE ==> qui convertiamo la dose in J/cm
    GafParameter(i,4)=i; %numero del GAF
    
    Dev_Standard_DO_GAF_i(i,1)=std(DO_GAF(:,i)); %dev sulla dose
  


end
  Dose_letta_GAF=GafParameter(:,1)./GafParameter(:,2); % cGy/pixel --> ? il valore medio di dose assorbita
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Calcolo l'errore sulla DENSITA' OTTICA dei GAF dello STACK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Derivata_parziale_y_stack=vpa(subs(dfy,y,StackParameter(:,5))); 
% Err_f_stack=sqrt((Derivata_parziale_x.^2).*((Dev_standard_Fondo+sensibility_intensity).^2)+(Derivata_parziale_y_stack.^2).*((StackParameter(:,6)+sensibility_intensity).^2)); 
% Err_DO_stack=2*StackParameter(:,7)./100; % calcolo l'errore del 2% sulla lettura dello scanner
% Err_tot_DO_Stack=Err_f_stack+Err_DO_stack;
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% Calcolo l'errore sulla DOSE dei GAF dello STACK usando la exp2 con termine noto 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% syms x a b c d e
% f= a*exp(x*b) + c*exp(x*d) + e;
% 
% dfx=diff(f,x); % derivata parziale rispetto a x
% dfa=diff(f,a); % derivata parziale rispetto a a
% dfb=diff(f,b); % derivata parziale rispetto a b
% dfc=diff(f,c); % derivata parziale rispetto a c
% dfd=diff(f,d); % derivata parziale rispetto a d
% dfe=diff(f,e); % derivata parziale rispetto a e
% 
% delta_a=1.407;
% delta_b=0.02583;
% delta_c=1.495;
% delta_d=0.01632;
% delta_e=1.511;
% 
% a1=575.7;
% b1=2.396;
% c1=-1557;
% d1=0.4782;
% 
% Derivata_parziale_e_Stack=vpa(subs(dfe,e,e));
% 
% syms f(x,c,d)
% f(x,c,d)=dfd;
% Derivata_parziale_d_Stack=vpa(subs(f(x,c1,d1),x,StackParameter(:,7))); % calcolo il valore numerico della derivata parziale rispetto a d per tutti i punti 
% 
% syms f(x,a,b)
% f(x,a,b)=dfb;
% Derivata_parziale_b_Stack=vpa(subs(f(x,a1,b1),x,StackParameter(:,7)));
% 
% syms f(x,d)
% f(x,d)=dfc;
% Derivata_parziale_c_Stack=vpa(subs(f(x,d1),x,StackParameter(:,7))); % calcolo il valore numerico della derivata parziale rispetto a c per tutti i punti x
% 
% syms f(x,b)
% f(x,b)=dfa;
% Derivata_parziale_a_Stack=vpa(subs(f(x,b1),x,StackParameter(:,7)));
% 
% syms f(x,a,b,c,d)
% f(x,a,b,c,d)=dfx;
% Derivata_parziale_x_Stack=vpa(subs(f(x,a,b,c,d),x,StackParameter(:,7))); % calcolo il valore numerico della derivata parziale rispetto a x per tutti i punti x
% 
% Err_DOSE_Stack=sqrt((Derivata_parziale_x_Stack.^2).*(Err_tot_DO_Stack.^2)+(Derivata_parziale_a_Stack.^2).*(delta_a.^2)+...
%     (Derivata_parziale_b_Stack.^2)*(delta_b.^2)+(Derivata_parziale_c_Stack.^2).*(delta_c.^2)+(Derivata_parziale_d_Stack.^2).*(delta_d.^2)+(Derivata_parziale_e_Stack.^2).*(delta_e.^2)); 
% 
% Err_DOSE_stack_J_cm(:,1)=Err_DOSE_Stack(:,1).*1.*0.00001*((PxcmFactor150).^2).*GafParameter(:,2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Conversione GAF ==> H2O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGAF=1:TotStack;
Spessore=NGAF.*(0.355);
%%
close all
figure(1)
hold on
plot(Spessore,Dose_letta_GAF,'LineWidth',1,'Color',[0 0 0]);
xlabel('Profondita equivalente di H2O [mm]');
ylabel('Dose [cGy/pixel]');
hold off

%%
%% TROVO  L'INDICE DEL GAF A CUI CORRISPONDE LA DOSE MASSIMA, QUELLA DA CUI PARTE LA CORREZIONE

massimoDose=max(Dose_letta_GAF);
[Max_index]=find(Dose_letta_GAF==massimoDose);
%%
%% RITAGLIO LE MATRICI OTTENUTE TRAMITE SIMULAZIONE SRIM

I=100-TotStack+1; %100 = il numero di righe e colonne di SRIM 

M=StoppingPowerMatrix.Factor41;
M1=M(I:100,I:100); %ritagliata

SM=StoppingPowerMatrix.StrimFactor;
SM1=SM(I:100,I:100);  %ritagliata

E=StoppingPowerMatrix.Energy; 
E1=E(1:TotStack,1);  %ritagliata

GafParameter(:,7)=E1(:,1); %l'energia dei protoni primari di SRIM


%%
%% NORMALIZZAZIONE E CORREZIONE ...????



%%==> Faccio la media della dose in tre GAF vicino l'enetrance e in cui
%%sono sicura che non c'? dipendenza dal LET

% ===>62MeV% 
Dose_Di_Normalizzazione=((GafParameter(6,1)+GafParameter(7,1)+GafParameter(8,1))/3)/GafParameter(1,2); % cGy/pixel
% ===>30MeV% 
%Dose_Di_Normalizzazione=GafParameter(1,1)/GafParameter(1,2);

% % % %Err_Dose_Di_Normalizzazione=sqrt((Err_DOSE_Stack(6,1)/3)^2+(Err_DOSE_Stack(7,1)/3)^2+(Err_DOSE_Stack(8,1)/3)^2);


% ==> Moltiplico tutte le righe della dose sperimentale per la dose media
% dell'entrance dello stack
DatiMarkus_Proporzionata(:,1)=DatiMarkus(:,3)*Dose_Di_Normalizzazione; % cGy/pixel

% % %Err_DatiMarkus_Proporzionata(:,1)=DatiMarkus(:,3)*Err_Dose_Di_Normalizzazione;

%% ==> Faccio una interpolazione dei dati sperimentali cosicch? da avere lo stesso numero di punti per lo stack e per la Markus+
Dose_Markus_interpol(:,1)=interp1(DatiMarkus(:,1), DatiMarkus_Proporzionata(:,1), Spessore,'linear');

% % % %Err_Dose_Markus_interpol=0;

%% divido la dose interpolata e moltiplicata per la dose media dello stack per la dose letta in tutto lo stack
%Fattore_di_correzione(:,1)=Dose_Markus_interpol./Dose_letta_GAF; % in cGy/pixel

Fattore_di_correzione(:,1)=Dose_Markus_interpol(:,1)./Dose_letta_GAF(:,1); % in cGy/pixel
for i=1:6
    Fattore_di_correzione(i,1)=1;
end

%%
%     a=Err_Dose_Markus_interpol(i,1)./Dose_letta_GAF_i(i,1);
% % %     a=Err_Dose_Markus_interpol./Dose_letta_GAF;
% % %     b=(Dose_Markus_interpol.*Err_DOSE_Stack)./((Dose_letta_GAF).^2);
% % %     Err_Fattore_di_correzione(:,1)=sqrt((a.^2)+(b.^2));

%%

%%==> In Gaf parameter corretto non c'? pi? la dose letta dal GAF ma una
%%dose opportunamente normalizzata per confrontarla con i dati sperimentali
GafParameter_Corretto=GafParameter(:,3).*Fattore_di_correzione;


%% CREO UN VETTORE COLONNA CHE CONTERRA' LA DOSE RILASCIATA SOLO DAI PROTONI PRIMARI DEL FASCIO DOPO LA CORREZIONE

L=TotStack-Max_index+1; %Maxindex=indice a cui corrisponde la dose massima => ci serve per caopire dov'? il picco 
%90-82+1 =9 perch? nel ciclo successivo mettiamo a zero gli ultimi 9 gaf
%dopo il picco

%carico le dosi dalla maggiore alla minore cio? dal picco all'entrance
%quindi inverto l'ordine con cui ? stata creata GafParameter

%Gaf=zeros(Tot,1); %creo un nuovo vettore
for i=L:TotStack %il ciclo va al contrario 
    
     DoseStack(i,1)=GafParameter((TotStack+1)-i,3); %metto in gaf quello che c'era in gafparameter ma con l'ordine inverso (partiamo da L e va fino Tot)
     %tot+1=90+1-(i=L=90-82=9)
     
     DoseStack_Corretta(i,1)=GafParameter_Corretto((TotStack+1)-i,1); %correggo per la sottostima dovuta al LET T 
end

%% SOTTRAZIONE PESATA

SottrazionePesata=zeros(TotStack,1); 
SottrazionePesata_Corretta=zeros(TotStack,1); %%questo ? lo stesso di prima ma la dose ? moltiplicata per il fattore di correzione

 for i=1:TotStack
       for j=1:TotStack
            SottrazionePesata(i,1)=SottrazionePesata(i,1)+(DoseStack(j,1)*M1(i,j)); %questa ? una sorta di sommatoria 
            SottrazionePesata_Corretta(i,1)=SottrazionePesata_Corretta(i,1)+(DoseStack_Corretta(j,1)*M1(i,j));
            
            
            %DoseStack(j,1)? il vettore invertito di DOSE che abbiamo creato
            %prima (convertita in J/cm cos? abbiamo la perdita di energia)
            
            %moltiplichiamo per M1 che ? il fattore di peso (cio? il
            %rapporto tra wk1/wi 

% S=StrimFactor
%la prima colonna ha in ogni riga la perdita di energia del protone che si
%ferma nel 90esimo gaf secondo quest'ordine: dal 90esimo gaf al primo
%quindi dal picco all'entrance

%M1 ? la matrice normalizzata ovvero ogni colonna di S ? divisa per il suo
%valore massimo

%j ? il contributo dei protoni secondari
            
       end
        
DoseStack(i,1)=(2.*DoseStack(i,1))-SottrazionePesata(i,1); %bisogna moltiplicare due volte gaf per far e tornare la formula
DoseStack_Corretta(i,1)=(2.*DoseStack_Corretta(i,1))-SottrazionePesata_Corretta(i,1);
 
 end

 %%
  

 GafParameter3_corretto(:,1)=GafParameter(:,3).*Fattore_di_correzione(:,1);
% % %  Err_GafParameter3_corretto(:,1)=((Err_DOSE_stack_J_cm(:,1).*Fattore_di_correzione(:,1)).^2+(GafParameter(:,3).*Err_Fattore_di_correzione(:,1)).^2).^(1/2);



%% Trasporto le dosi corrette nella conolla 5 di GafParameter e sistemo l'ordine

%%==> Riempio la quinta colonna della matrice GafParameter

GafParameter5_corretto(:,1)=zeros(TotStack,1); % con correzione della sottostima 
% % ERR_GafParameter5(:,1)=Err_DOSE_stack_J_cm(:,1);
% % ERR_GafParameter5_corretto(:,1)=Err_GafParameter3_corretto(:,1);

for i=1:Max_index
    GafParameter(i,5)=DoseStack((TotStack+1)-i); 
    GafParameter5_corretto(i,1)=DoseStack_Corretta((TotStack+1)-i); 
end

for i=Max_index+1:TotStack
    GafParameter(i,5)=GafParameter(i,3);
    GafParameter5_corretto(i,1)=GafParameter_Corretto(i,1); 
end





%% %% PLOTTO LA DISTRIBUZIONE DI ENERGIA IN PROFONDITA' DOVUTA SOLO AI PROTONI PRIMARI
figure(20)
hold on
plot(Spessore, GafParameter(:,5), 'r-');
%plot(Spessore, GafParameter(:,3), 'r-');
%plot(Spessore, GafParameter(:,1), 'r-');
%errorbar(Spessore, GafParameter(:,5),ERR_GafParameter5(:,1),'go');
plot(Spessore, GafParameter5_corretto, 'k-');
%errorbar(Spessore, GafParameter5_corretto(:,1),ERR_GafParameter5_corretto(:,1),'bo');
legend('energia rilasciata dai soli protoni primari','energia rilasciata dai soli protoni primari con correzione della sottostima');
title('Confronto della distribuzione di energia in profondit? dovuta solo ai protoni primari di un fascio di 60 MeV FE: con correzione della sottostima e senza');
xlabel('mmH2O');
ylabel('Stopping Power [J/cm]');
xlim([0 32]);
hold off


%% TROVO IL NUMERO DEI PROTONI CHE SI SONO FERMATI IN OGNI GAF ALL'INTERNO DELLA ROI

% % % % ERR_GafParameter6(:,1)=zeros(TotStack,1);
GafParameter6_corretto(:,1)=zeros(TotStack,1); % con correzione della sottostima 
% % % % ERR_GafParameter6_corretto(:,1)=zeros(TotStack,1);
%%
for j=1:TotStack
      i=(TotStack+1)-j;
      
      GafParameter(j,6)=GafParameter(j,5)./(SM1(i,i)*1.6022*(10^-11)); % trovo #protoni e converto la matrice SRIM in J/cm
      %ERR_GafParameter6(j,1)=ERR_GafParameter5(j,1)./(SM1(i,i)*1.6022*(10^-11));  
      GafParameter6_corretto(j,1)=GafParameter5_corretto(j,1)./(SM1(i,i)*1.6022*(10^-11));  % trovo #protoni e converto la matrice SRIM in J/cm
% % %       ERR_GafParameter6_corretto(j,1)=ERR_GafParameter5_corretto(j,1)/(SM1(i,i)*1.6022*(10^-11));  
      % ===> applico la formula della dose per ottenere 
      
      
end
%%


         Fluenza_cm2=GafParameter(:,6)./AreaRectStack;
         
% %          a=ERR_GafParameter6./AreaRectStack;
% %          b=(GafParameter(:,6).*Err_AreaRectStack)./(AreaRectStack).^2;
% %          ERR_Fluenza_cm2=sqrt((a.^2)+(b.^2));
% %          
         Fluenza_cm2_corretta= GafParameter6_corretto./AreaRectStack;
         
% % %          c=ERR_GafParameter6_corretto./AreaRectStack;
% % %          d=(GafParameter6_corretto.*Err_AreaRectStack)./(AreaRectStack).^2;
% % %          ERR_Fluenza_cm2_corretta=sqrt((c.^2)+(d.^2));
% % %          
         Fluenza_px(:,1)=GafParameter(:,6)/GafParameter(1,1);
         Fluenza_px_corretta(:,1)=GafParameter6_corretto(:,1)/GafParameter(1,1);
         
% % %          e=ERR_GafParameter6_corretto./GafParameter(1,1);
% % %          f=(GafParameter6_corretto.*1)./(GafParameter(1,1)).^2;
% % %          ERR_Fluenza_px_corretta=sqrt((e.^2)+(f.^2));
% % %       

%%

% % NGAF=1:Tot;
% % Spessore=NGAF.*(0.355);

xmax= 60;

close all
figure(1)
subplot(3,3,1)
plot(Spessore,GafParameter(:,1));
ylabel('Gy')
xlim([0 xmax]);
title('Somma delle dosi in tutti i pixel della ROI');

subplot(3,3,2)
plot(Spessore,GafParameter(:,2));
ylabel('a.u.')
xlim([0 xmax]);
title('Numero di pixel della ROI');

subplot(3,3,3)
plot(Spessore,GafParameter(:,3));
ylabel('J/cm')
xlim([0 xmax]);
title('Somma delle dosi in tutti i pixel della ROI');

subplot(3,3,4)
plot(Spessore,GafParameter(:,4));
ylabel('a.u.')
xlim([0 xmax]);
title('Numero del Gaf');

subplot(3,3,5)
plot(Spessore,GafParameter(:,5),Spessore, GafParameter5_corretto);
ylabel('J/cm')
xlim([0 xmax]);
title('Distribuzione energetica dei protoni incidenti');


subplot(3,3,6)
plot(Spessore,GafParameter(:,6),Spessore,GafParameter6_corretto);
ylabel('a.u.')
xlim([0 xmax]);
title('Numero di protoni che si sono fermati in ogni Gaf');

subplot(3,3,7)
plot(Spessore,GafParameter(:,7));
ylabel('MeV')
xlim([0 xmax]);
title('Energia dei protoni che si sono fermati in ogni Gaf');

subplot(3,3,8)
plot(Spessore,Dose_letta_GAF,'LineWidth',1,'Color',[0 0 0]);
ylabel('cGy/pixel');
xlim([0 xmax]);
title('Dose in ogni pixel');

subplot(3,3,9)
plot(Spessore,Fluenza_cm2,'LineWidth',1,'Color',[0 0 0]);
ylabel('a.u.');
xlim([0 xmax]);
title('Fluenza al cm2');

figure
plot(GafParameter(:,7),GafParameter(:,5))


