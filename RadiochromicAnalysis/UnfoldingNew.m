% % % 1 - Procedura di Unfolding: il codice originale deriva da una versione 
% % %     non funzionante di Pietro Pisciotta;
% % % 2 - La versione funzionante dell'intero codice di riferimento
% % %     (originariamente calibrazione + unfolding + calcolo e propagazione degli errori su tutto) 
% % %     è stato il lavoro di Tesi di Cristina - anno 2019
% % % 3 - Nell'estate 2021 Giada rende più compatto il codice
% % % 4 - Correzione del codice sul calcolo della dose (Cristina) e versione attuale senza calibrazione - 2022

%% 
clear all
close all

%% Load of the stopping power matrix and definition of some constant values

MatrixSP=load('40umAl_HDV2_70MeV.mat'); 

S='HDV20';
Type='.tif';
Zero='0';

DPI=600;
PxcmFactor=2.54/DPI; %dimensione di un pixel in cm
Spot_cm2=pi*(0.75^2); % sezione trasversale del fascio in cm2 usando un collimatore con diametro da 15 mm
Spot_px=Spot_cm2/(PxcmFactor^2); % sezione trasversale del fascio in #pixel usando un collimatore con diametro da 15 mm
sensibility_intensity=1; % è la sensibilità di misura con cui valuto i livelli di intensit?

%% imposto una curva definendo i parametri per il fit

a=67958.206583;
b=4856.320357;
c=6906.140299;

poly3_0=@(x) a*x.^(3)+b*x.^(2)+c*x;
% [fitobject,gof] = fit(OD(:,1),DOSE(:,1),'poly3_0'); 
PR=poly3_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calcolo intensità e dose media in ogni ROI di ogni gaf dello stack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Seleziona una immagine qualsiasi nella cartella dei dati che vuoi analizzare')
[FileName_Stack,PathName_Stack]=uigetfile('*.tif','Seleziona una immagine qualsiasi nella cartella dei dati che vuoi analizzare.');
TotStack=input('\n Inserisci il numero di gaf dello stack\n');
fullchosenfileSTACK = [PathName_Stack FileName_Stack];

figure(1) 
RCF=imread(fullchosenfileSTACK);
[R, rect_stack]=imcrop(RCF);

Numero_di_pixel_ROI_Stack=round(rect_stack(1,3))*round(rect_stack(1,4));
base=round(rect_stack(1,3))*PxcmFactor; % cm
altezza=round(rect_stack(1,4))*PxcmFactor; % cm
AreaRectStack=base*altezza; % area della ROI in cm2 
Err_AreaRectStack=sqrt((base)^2+(altezza)^2)*PxcmFactor; % errore della misura dell'area della ROI in cm2

disp('Seleziona immagine da usare come fondo per lo stack')
[FileName_Fondo,PathName_Fondo]=uigetfile('*.tif','Seleziona immagine da usare come fondo');
fullchosenfile_Fondo = [PathName_Fondo FileName_Fondo];
Fondo_Stack=imread(fullchosenfile_Fondo);

figure(2)
[R rect_fondo_stack]=imcrop(Fondo_Stack);
AreaRect_Fondo_Stack=[round(rect_fondo_stack(1,3))*PxcmFactor*round(rect_fondo_stack(1,4))*PxcmFactor]; %area della ROI in cm2 

Numero_di_pixel_ROI_Fondo_Stack=round(rect_fondo_stack(1,3))*round(rect_fondo_stack(1,4));

%% Lavoro sul fondo

Vettore_Intensita_Fondo_Stack=zeros(Numero_di_pixel_ROI_Fondo_Stack,1);
t=1;
for j=round(rect_fondo_stack(1,1)):round(rect_fondo_stack(1,1))+round(rect_fondo_stack(1,3)) 
    for k=round(rect_fondo_stack(1,2)):round(rect_fondo_stack(1,2))+round(rect_fondo_stack(1,4))
        Vettore_Intensita_Fondo_Stack(t,1)=Fondo_Stack(round(k),round(j),1);
        t=t+1;
    end
end

Dev_standard_Fondo_Stack=std(Vettore_Intensita_Fondo_Stack);
Average_Fondo_Stack=mean(Vettore_Intensita_Fondo_Stack);

%% Lavoro sullo stack

StackParameter=zeros(TotStack,8);
GafParameter=zeros(TotStack,7);
DO_GAF=zeros(Numero_di_pixel_ROI_Stack,TotStack);
Dev_Standard_DO_GAF=zeros(TotStack,1);

for i=1:TotStack 
    if i<10    
       AT=int2str(i);
       A=strcat(Zero,AT);
    else 
       A=int2str(i);
    end

    FileName_Stack=strcat(S,A,Type);
    fullchosenfileStack = [PathName_Stack FileName_Stack];
    GAFSTACK=imread(fullchosenfileStack);
   
    figure(3)
    imshow(GAFSTACK);

    Xmin=round(rect_stack(1,1));
    Xmax=round(rect_stack(1,1)+rect_stack(1,3));
    Ymin=round(rect_stack(1,2));
    Ymax=round(rect_stack(1,2)+rect_stack(1,4));
    
    t=1;
    m=1;

    for j=Xmin:Xmax
        for k=Ymin:Ymax
            PVR=double(GAFSTACK(k,j,1));
                   
            Livello_Minimo=Average_Fondo_Stack;
            
            if PVR<Livello_Minimo % check per capire se ci sono macchie 
               
               Livello_Minimo=PVR;
               Matrice_Intensita_Gaf_Stack(t,i)=PVR;
               StackParameter(i,1)=StackParameter(i,1)+PVR; 
               StackParameter(i,2)=StackParameter(i,2)+1;
               t=t+1;
               DO_PVR=log10(Average_Fondo_Stack/PVR);
               DO_GAF(m,i)=DO_PVR;    
               DoseSinglePixel= PR(DO_PVR); 
               GafParameter(i,1)=GafParameter(i,1)+DoseSinglePixel; %qui abbiamo la somma della dose di ogni pixel della roi
               GafParameter(i,2)=GafParameter(i,2)+1; %qui abbiamo il numero di pixel
               m=m+1;
            
            else
                
               Matrice_Intensita_Gaf_Stack(t,i)=PVR; 
               StackParameter(i,1)=StackParameter(i,1)+PVR; %PVR per ogni ROI di ogni GAF
               StackParameter(i,2)=StackParameter(i,2)+1; %numero di pixel della ROI

               t=t+1;
               DO_PVR=log10(Average_Fondo_Stack/PVR); %qui calcoliamo la densita ottica del singolo pixel mentre in StackParameter(i,7) calcoliamo il valor medio
               DO_GAF(m,i)=DO_PVR;    
               DoseSinglePixel= PR(DO_PVR); 
               GafParameter(i,1)=GafParameter(i,1)+DoseSinglePixel; %qui abbiamo la somma della dose di ogni pixel della roi
               GafParameter(i,2)=GafParameter(i,2)+1; %qui abbiamo il numero di pixel
               m=m+1;
            end
        end
    end
    
    StackParameter(i,3)=Livello_Minimo; 
    StackParameter(i,4)=i; %numero del GAF
    StackParameter(i,5)=StackParameter(i,1)./StackParameter(i,2); %  valore medio dei PVR della ROI
    StackParameter(i,6)=std(Matrice_Intensita_Gaf_Stack(:,i))./sqrt(Numero_di_pixel_ROI_Stack); 
    StackParameter(i,7)=log10(Average_Fondo_Stack/StackParameter(i,5)); %  DO media, ovvero calcolata sulla media dell'intensita della ROI 
    StackParameter(i,8)=PR(StackParameter(i,7)); % la DOSE media in cGy calcolata considerando la DO media 
    
    
    GafParameter(i,3)=(GafParameter(i,1)*1*0.00001)*((PxcmFactor)^2);  %% converto la dose in J/cm
    GafParameter(i,4)=i; % numero del GAF
    Dev_Standard_DO_GAF_i(i,1)=std(DO_GAF(:,i)); % dev sulla dose
  
end
  
Dose_letta_GAF=GafParameter(:,1)./GafParameter(:,2); % cGy/pixel --> è il valore medio di dose assorbita

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calcolo l'errore sulla DENSITA' OTTICA dei GAF dello STACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derivata_parziale_y_stack=vpa(subs(dfy,y,StackParameter(:,5))); 
% Err_f_stack=sqrt((Derivata_parziale_x.^2).*((Dev_standard_Fondo_Stack+sensibility_intensity).^2)+(Derivata_parziale_y_stack.^2).*((StackParameter(:,6)+sensibility_intensity).^2)); 
% Err_DO_stack=2*StackParameter(:,7)./100; % calcolo l'errore del 2% sulla lettura dello scanner
% Err_tot_DO_Stack=Err_f_stack+Err_DO_stack;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calcolo l'errore sulla DOSE dei GAF dello STACK usando la exp2 con termine noto 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Err_DOSE_stack_J_cm(:,1)=Err_DOSE_Stack(:,1).*1.*0.00001*((PxcmFactor).^2).*GafParameter(:,2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Conversione GAF ==> H2O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Strato_Attivo=195; % um di acqua
% s=0;
% for i=1:TotStack
%     Spessore(i,1)=Strato_Attivo+s*355;  % per EBT3
%     s=s+1;
% end

Strato_Attivo=15; % um di acqua
s=0;
for i=1:TotStack
    Spessore(i,1)=Strato_Attivo+s*150;    % per HDV2
    s=s+1;
end

%%
close all

figure(1)
hold on
% plot(Spessore,Dose_letta_GAF./Dose_letta_GAF(1,1),Markus_FE_21_02_2019(:,1),Markus_FE_21_02_2019(:,2)./(Markus_FE_21_02_2019(1,2)),'LineWidth',1,'Color',[0 0 0]);
% plot(Spessore,Dose_letta_GAF./Dose_letta_GAF(1,1),'LineWidth',1,'Color',[0 0 0]);
plot(Spessore,Dose_letta_GAF,'bo');
title('Dose distribution of the primary protons in water');
xlabel('Profondità equivalente di H2O [um]');
ylabel('Dose [cGy/pixel]');
hold off

%% TROVO  L'INDICE DEL GAF A CUI CORRISPONDE LA DOSE MASSIMA, QUELLA DA CUI PARTE LA CORREZIONE

massimoDose=max(Dose_letta_GAF);
[Max_index]=find(Dose_letta_GAF==massimoDose);

%% RITAGLIO LE MATRICI OTTENUTE TRAMITE SIMULAZIONE SRIM

N=size(MatrixSP.Matrix_Normalized,1);
I=N-TotStack+1;  

M=MatrixSP.Matrix_Normalized;
M1=M(I:N,I:N); %ritagliata

SM=MatrixSP.Matrix;
SM1=SM(I:N,I:N);  %ritagliata

E=MatrixSP.Energia_MeV; 
E1=E(1:TotStack,1);  %ritagliata
GafParameter(:,7)=E1(:,1); 

%% NORMALIZZAZIONE E CORREZIONE ...????

%%==> Faccio la media della dose in tre GAF vicino l'enetrance e in cui
%%sono sicura che non c'? dipendenza dal LET
% Dose_Di_Normalizzazione=((GafParameter(6,1)+GafParameter(7,1)+GafParameter(8,1))/3)/GafParameter(1,2); % cGy/pixel
% Err_Dose_Di_Normalizzazione=sqrt((Err_DOSE_Stack(6,1)/3)^2+(Err_DOSE_Stack(7,1)/3)^2+(Err_DOSE_Stack(8,1)/3)^2);


% ==> Moltiplico tutte le righe della dose sperimentale per la dose media
% dell'entrance dello stack
% Markus_FE_21_02_2019_Proporzionata(:,1)=Markus_FE_21_02_2019(:,3)*Dose_Di_Normalizzazione; % cGy/pixel
% Err_Markus_FE_21_02_2019_Proporzionata(:,1)=Markus_FE_21_02_2019(:,3)*Err_Dose_Di_Normalizzazione;

%% ==> Faccio una interpolazione dei dati sperimentali cosicch? da avere lo stesso numero di punti per lo stack e per la Markus+
% Dose_Markus_interpol(:,1)=interp1(Markus_FE_21_02_2019(:,1), Markus_FE_21_02_2019_Proporzionata(:,1), Spessore,'linear');
% Err_Dose_Markus_interpol=0;

%% divido la dose interpolata e moltiplicata per la dose media dello stack per la dose letta in tutto lo stack
%Fattore_di_correzione(:,1)=Dose_Markus_interpol./Dose_letta_GAF; % in cGy/pixel

% Fattore_di_correzione(:,1)=Dose_Markus_interpol(:,1)./Dose_letta_GAF(:,1); % in cGy/pixel
% for i=1:6
%     Fattore_di_correzione(i,1)=1;
% end

%%
%     a=Err_Dose_Markus_interpol(i,1)./Dose_letta_GAF_i(i,1);
%     a=Err_Dose_Markus_interpol./Dose_letta_GAF;
%     b=(Dose_Markus_interpol.*Err_DOSE_Stack)./((Dose_letta_GAF).^2);
%     Err_Fattore_di_correzione(:,1)=sqrt((a.^2)+(b.^2));

%%

%%==> In Gaf parameter corretto non c'? pi? la dose letta dal GAF ma una
%%dose opportunamente normalizzata per confrontarla con i dati sperimentali

% GafParameter_Corretto=GafParameter(:,3).*Fattore_di_correzione;


%% CREO UN VETTORE COLONNA CHE CONTERRA' LA DOSE RILASCIATA SOLO DAI PROTONI PRIMARI DEL FASCIO DOPO LA CORREZIONE

L=TotStack-Max_index+1; 

for i=L:TotStack %il ciclo va al contrario 
    
     DoseStack(i,1)=GafParameter((TotStack+1)-i,3); % metto in gaf quello che c'era in gafparameter ma con l'ordine inverso (partiamo da L e va fino Tot)
%      DoseStack_Corretta(i,1)=GafParameter_Corretto((TotStack+1)-i,1); % correggo per la sottostima dovuta al LET T 
end

%% SOTTRAZIONE PESATA

SottrazionePesata=zeros(TotStack,1); 
% SottrazionePesata_Corretta=zeros(TotStack,1); 

 for i=1:TotStack
       for j=1:TotStack
            SottrazionePesata(i,1)=SottrazionePesata(i,1)+(DoseStack(j,1)*M1(i,j)); 
%             SottrazionePesata_Corretta(i,1)=SottrazionePesata_Corretta(i,1)+(DoseStack_Corretta(j,1)*M1(i,j));
            
            
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
        
DoseStack(i,1)=(2.*DoseStack(i,1))-SottrazionePesata(i,1); % bisogna moltiplicare due volte gaf per far e tornare la formula
% DoseStack_Corretta(i,1)=(2.*DoseStack_Corretta(i,1))-SottrazionePesata_Corretta(i,1);
 
 end
  
%  GafParameter3_corretto(:,1)=GafParameter(:,3).*Fattore_di_correzione(:,1);
%  Err_GafParameter3_corretto(:,1)=((Err_DOSE_stack_J_cm(:,1).*Fattore_di_correzione(:,1)).^2+(GafParameter(:,3).*Err_Fattore_di_correzione(:,1)).^2).^(1/2);

%% Trasporto le dosi corrette nella conolla 5 di GafParameter e sistemo l'ordine

% GafParameter5_corretto(:,1)=zeros(TotStack,1); % con correzione della sottostima 
% ERR_GafParameter5(:,1)=Err_DOSE_stack_J_cm(:,1);
% ERR_GafParameter5_corretto(:,1)=Err_GafParameter3_corretto(:,1);

for i=1:Max_index
    GafParameter(i,5)=DoseStack((TotStack+1)-i); 
%     GafParameter5_corretto(i,1)=DoseStack_Corretta((TotStack+1)-i); 
end

for i=Max_index+1:TotStack
    GafParameter(i,5)=GafParameter(i,3);
%     GafParameter5_corretto(i,1)=GafParameter_Corretto(i,1); 
end


%% %% PLOTTO LA DISTRIBUZIONE DI ENERGIA IN PROFONDITA' DOVUTA SOLO AI PROTONI PRIMARI
% figure(2)
% hold on
% plot(Spessore, GafParameter(:,5), 'r-');
% %plot(Spessore, GafParameter(:,3), 'r-');
% %plot(Spessore, GafParameter(:,1), 'r-');
% %errorbar(Spessore, GafParameter(:,5),ERR_GafParameter5(:,1),'go');
% % plot(Spessore, GafParameter5_corretto, 'k-');
% %errorbar(Spessore, GafParameter5_corretto(:,1),ERR_GafParameter5_corretto(:,1),'bo');
% % legend('energia rilasciata dai soli protoni primari','energia rilasciata dai soli protoni primari con correzione della sottostima');
% title('Confronto della distribuzione di energia in profondit? dovuta solo ai protoni primari di un fascio di 60 MeV FE: con correzione della sottostima e senza');
% xlabel('mmH2O');
% ylabel('Stopping Power [J/cm]');
% % xlim([0 32]);
% hold off

%% TROVO IL NUMERO DEI PROTONI CHE SI SONO FERMATI IN OGNI GAF ALL'INTERNO DELLA ROI

% ERR_GafParameter6(:,1)=zeros(TotStack,1);
% GafParameter6_corretto(:,1)=zeros(TotStack,1); % con correzione della sottostima 
% ERR_GafParameter6_corretto(:,1)=zeros(TotStack,1);

for j=1:TotStack
      i=(TotStack+1)-j;
      
      GafParameter(j,6)=GafParameter(j,5)./(SM1(i,i)*1.6022*(10^-11)); % trovo #protoni e converto la matrice SRIM in J/cm
      %ERR_GafParameter6(j,1)=ERR_GafParameter5(j,1)./(SM1(i,i)*1.6022*(10^-11));  
%       GafParameter6_corretto(j,1)=GafParameter5_corretto(j,1)./(SM1(i,i)*1.6022*(10^-11));  % trovo #protoni e converto la matrice SRIM in J/cm
%       ERR_GafParameter6_corretto(j,1)=ERR_GafParameter5_corretto(j,1)/(SM1(i,i)*1.6022*(10^-11));  
      % ===> applico la formula della dose per ottenere      
end

%%
         Fluenza_cm2=GafParameter(:,6)./AreaRectStack;
         
%          a=ERR_GafParameter6./AreaRectStack;
%          b=(GafParameter(:,6).*Err_AreaRectStack)./(AreaRectStack).^2;
%          ERR_Fluenza_cm2=sqrt((a.^2)+(b.^2));
         
%          Fluenza_cm2_corretta= GafParameter6_corretto./AreaRectStack;
%          
%          c=ERR_GafParameter6_corretto./AreaRectStack;
%          d=(GafParameter6_corretto.*Err_AreaRectStack)./(AreaRectStack).^2;
%          ERR_Fluenza_cm2_corretta=sqrt((c.^2)+(d.^2));
         
         Fluenza_px(j,1)=GafParameter(j,6)/GafParameter(1,1);
%          Fluenza_px_corretta(j,1)=GafParameter6_corretto(j,1)/GafParameter(1,1);
         
%          e=ERR_GafParameter6_corretto./GafParameter(1,1);
%          f=(GafParameter6_corretto.*1)./(GafParameter(1,1)).^2;
%          ERR_Fluenza_px_corretta=sqrt((e.^2)+(f.^2));

%% PLOTTO LA FLUENZA

figure(9)
hold on
plot(GafParameter(:,7),Fluenza_cm2(:,1),'ko');
xlabel('Energy [MeV]');
ylabel('Fluence [protons/cm^2]');
title('Fluence from deconvolution');
hold off
