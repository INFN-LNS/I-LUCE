%%% Mariacristina Guarrera 2025

clear all
close all

%% Add folders

addpath("Functions/");

%% Load of the stopping power matrix and definition of the film name

MatrixSP=load('0umAl_EBT3_70MeV.mat'); 

disp("Enter the name of the RCF type:  ");
disp('For example: EBT3, EBT4, HDV2');
disp('NB: Valid for full-GAF stack, i.e. without layers in between or a mixture of GAF types');
S = input('Enter answer: ', 's');
Type='.tif';
Zero='0';
ZeroZero='00';

%% Asks for the DPI value and definition of the pixel dimension in cm

prompt = 'Enter the DPI value used to scan the RCF films: '; % DPI (dots per inch) of the scanned image
DPI = input(prompt); % usually, it's 150 for LNS films and 600 for ELI films
PxcmFactor=2.54/DPI; % pixel dimension in cm

%% Upload a valid Calibration

disp('Have you a calibration file already prepared?')
disp('y/n or press x to exit')

P = input('Enter answer: ', 's');

    switch P
        case 'y'
            disp('Select the file with the calibration data and fit parameters.')
            filter = {'*.txt'};
            [FileNameCalibrationDataAndFit, path] = uigetfile(filter,'File Selection','CalibrationDataAndFit.txt');
            pathForCalibrationDataAndFit = strcat(path, FileNameCalibrationDataAndFit);
            [P1, P2, P3, BackgroundValue, Doses, NetOpticalDensities] = ReadCalibrationDataAndFit(pathForCalibrationDataAndFit);
       
        case 'n'
            disp('Enter the three fit parameters in the order P1, P2 and P3')
            FitParameterAsInput = [];
            for i = 1:1:3
                FitParameter = input('Insert the fit parameter P' + string(i) + ':');
                FitParametersAsInput(i) = FitParameter;
            end

            P1=FitParametersAsInput(1,1); 
            P2=FitParametersAsInput(1,2);  
            P3=FitParametersAsInput(1,3);  

        case 'x'
            disp('The program has been terminated')
        otherwise
            disp('Insert a valid value')
            return
    end
    
poly3_0=@(x) P3*x.^(3)+P2*x.^(2)+P1*x;
PR=poly3_0;

%% Seleziono lo stack e definisco una ROI

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

%% Scelgo l'immagine da usare come sfondo per lo stack e definisco una ROI

disp('Seleziona immagine da usare come fondo per lo stack')
[FileName_Fondo,PathName_Fondo]=uigetfile('*.tif','Seleziona immagine da usare come fondo');
fullchosenfile_Fondo = [PathName_Fondo FileName_Fondo];
Fondo_Stack=imread(fullchosenfile_Fondo);

figure(2)
[R rect_fondo_stack]=imcrop(Fondo_Stack);
AreaRect_Fondo_Stack=[round(rect_fondo_stack(1,3))*PxcmFactor*round(rect_fondo_stack(1,4))*PxcmFactor]; %area della ROI in cm2 
Numero_di_pixel_ROI_Fondo_Stack=round(rect_fondo_stack(1,3))*round(rect_fondo_stack(1,4));

%% Calcolo valore medio e dev standard dell'intensità dei pixel dentro la ROI del fondo dello stack

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

GafParameter=zeros(TotStack,8);
Matrice_Intensita_Gaf_Stack=zeros(Numero_di_pixel_ROI_Stack,TotStack);
Matrice_DO_GAF=zeros(Numero_di_pixel_ROI_Stack,TotStack);
Dev_Standard_DO_GAF=zeros(TotStack,1);

for i=1:TotStack 

    if      i<10    
            AT=int2str(i);
            A=strcat(ZeroZero,AT);
    elseif  (i>=10)&&(i<100)
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

    for j=Xmin:Xmax % colonne
        for k=Ymin:Ymax % righe
            PVR=double(GAFSTACK(k,j,1)); % intensità nel rosso del singolo pixel
               
               Matrice_Intensita_Gaf_Stack(t,i)=PVR; 
               t=t+1;
  
               Matrice_DO_GAF(m,i)=log10(Average_Fondo_Stack/PVR); % densita ottica del singolo pixel
               GafParameter(i,1)=GafParameter(i,1)+PR(Matrice_DO_GAF(m,i)); % qui abbiamo la somma della dose di ogni pixel della roi
               GafParameter(i,2)=GafParameter(i,2)+1; % qui abbiamo il numero di pixel
               GafParameter(i,8)=GafParameter(i,8)+PVR; % somma dei PVR per ogni ROI di ogni GAF
               m=m+1;
        end
    end
       
    GafParameter(i,3)=(GafParameter(i,1)*1*0.00001)*((PxcmFactor)^2);  %% converto la dose in J/cm
    GafParameter(i,4)=i; % numero del GAF
    Dev_Standard_DO_GAF_i(i,1)=std(Matrice_DO_GAF(:,i)); % dev sulla dose
  
end
  
Dose_letta_GAF=GafParameter(:,1)./GafParameter(:,2); % Gy/pixel --> è il valore medio di dose assorbita

%% creation of depth-vector in um of water

switch S
    case "EBT3"
    
        Strato_Attivo=0.195; % mm 
        s=0;
        for i=1:TotStack
            Spessore(i,1)=Strato_Attivo+s*0.355;
            s=s+1;
        end

    case "EBT4"
    
        Strato_Attivo=0.195; % mm 
        s=0;
        for i=1:TotStack
            Spessore(i,1)=Strato_Attivo+s*0.355;
            s=s+1;
        end

    case "HDV2"

        Strato_Attivo=0.015; % mm 
        s=0;
        for i=1:TotStack
            Spessore(i,1)=Strato_Attivo+s*0.150;   
            s=s+1;
        end
    otherwise 
        disp('Insert a valid input')
        return
end 

%%

figure(10)
hold on
plot(Spessore,Dose_letta_GAF,'b.');
title('Dose distribution of the protons in water');
xlabel('Profondità equivalente di H2O [mm]');
ylabel('Dose [Gy/pixel]');
hold off

%% TROVO  L'INDICE DEL GAF A CUI CORRISPONDE LA DOSE MASSIMA, QUELLA DA CUI PARTE LA CORREZIONE

massimoDose=max(Dose_letta_GAF);
[Max_index]=find(Dose_letta_GAF==massimoDose);
Max_index

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

%% CREO UN VETTORE COLONNA CHE CONTERRA' LA DOSE RILASCIATA SOLO DAI PROTONI PRIMARI DEL FASCIO 

L=TotStack-Max_index+1; 

for i=L:TotStack % il ciclo va al contrario 
    DoseStack(i,1)=GafParameter((TotStack+1)-i,3); % metto in DoseStack quello che c'era in gafparameter (:,3) ma con l'ordine inverso (partiamo da L, fino a TotStack)
end

%% SOTTRAZIONE PESATA

SottrazionePesata=zeros(TotStack,1); 

for i=1:TotStack
       for j=1:TotStack
            SottrazionePesata(i,1)=SottrazionePesata(i,1)+(DoseStack(j,1)*M1(i,j));      
       end
       DoseStack(i,1)=(2.*DoseStack(i,1))-SottrazionePesata(i,1); % bisogna moltiplicare due volte gaf per far e tornare la formula
end

%% Trasporto le dosi corrette nella colonna 5 di GafParameter e sistemo l'ordine

for i=1:Max_index
    GafParameter(i,5)=DoseStack((TotStack+1)-i);  
end

for i=Max_index+1:TotStack
    GafParameter(i,5)=GafParameter(i,3); 
end


%% %% PLOTTO LA DISTRIBUZIONE DI ENERGIA IN PROFONDITA' DOVUTA SOLO AI PROTONI PRIMARI

figure(20)
hold on
plot(Spessore, GafParameter(:,5), 'r-');
title('Distribuzione di energia dei soli protoni primari');
xlabel('mm H2O');
ylabel('Stopping Power [J/cm]');
xlim([0 35]);
hold off

%% TROVO IL NUMERO DEI PROTONI CHE SI SONO FERMATI IN OGNI GAF ALL'INTERNO DELLA ROI

for j=1:TotStack
      i=(TotStack+1)-j;
      GafParameter(j,6)=GafParameter(j,5)./(SM1(i,i)*1.6022*(10^-11)); % trovo #protoni e converto la matrice SRIM in J/cm     
end

%% Calcolo e plotto la fluenza 

Fluenza_cm2=GafParameter(:,6)./AreaRectStack;
Fluenza_px(j,1)=GafParameter(j,6)/GafParameter(1,1);

%%

figure(30)
hold on
plot(GafParameter(:,7),Fluenza_cm2(:,1),'k.');
xlabel('Energy [MeV]');
ylabel('Fluence [protons/cm^2]');
title('Fluence from deconvolution');
hold off
