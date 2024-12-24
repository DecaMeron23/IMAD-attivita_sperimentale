%% Componenti del Gruppo

% Matricola: 1053435 --- Nome: Raffaele Giacomo Giovanni Di Maio
% Matricola: 1080976 --- Nome: Emilio Meroni

%% Inizzializzazione
clc
clear
close all
rng('default');

set(0, 'defaultAxesFontSize', 12)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultlegendInterpreter','latex')

addpath("Dati\");

%% import data
dataAll = load("datamotor_2024-12-06_16_53_33.mat");
t_raw_step = dataAll.t; % time [s]
u_raw_step = dataAll.u; % input:  PWM duty cycle [0-100%]
y_raw_step = dataAll.y; % output: motor speed [rpm]

dataAll = load("datamotor_2024-12-06_17_18_14.mat");
t_raw_noise = dataAll.t; % time [s]
u_raw_noise = dataAll.u; % input:  PWM duty cycle [0-100%]
y_raw_noise = dataAll.y; % output: motor speed [rpm]

%% Studio preliminare
Ts = t_raw_step(2) - t_raw_step(1);

plot(t_raw_step , y_raw_step);
hold on
plot(t_raw_step , u_raw_step);

xlim([11.0640 11.1290])
ylim([0.0 54.9])

ax = gca;
chart = ax.Children(1);
datatip(chart,11.08,50);
chart = ax.Children(2);
datatip(chart,11.11,20.23);

hold off

figure()

plot(t_raw_noise, y_raw_noise);
hold on
plot(t_raw_noise , u_raw_noise);
hold off

%si osserva un ritardo puro di 5 passi
ritardo = 5;

%% Pulizia segnali

intervallo = floor(6/Ts) : floor(23/Ts);

t_noise = t_raw_noise(intervallo); % time [s]
u_noise = u_raw_noise(intervallo); % input:  PWM duty cycle [0-100%]
y_noise = y_raw_noise(intervallo);

u_media = mean(u_noise);
u_noise = u_noise - u_media;

y_media = mean(y_noise);
y_noise = y_noise - y_media;

t_noise = t_noise - t_noise(1); 

figure
plot(t_noise , u_noise);
hold on
plot(t_noise , y_noise);

%% Identificazione non parametrica

data = iddata(y_noise, u_noise, Ts);
FDT_non_parametrica = spa(data);

figure
h = bodeplot(FDT_non_parametrica);
showConfidence(h,3)

%% Identificazione parametrica
% dal diagramma di bode del modulo ipotizziamo una pulsazione critica pari
% a 30 rad/s che corrisponde a circa 5hz quindi ricampioniamo il segnale
% a circa 100hz (10 w_c/2pi < fs < 30 w_c/2pi)

% riduciamo di un fattore 2
t_noise_dec = decimate(t_noise , 2);
y_noise_dec = decimate(y_noise , 2);
u_noise_dec = decimate(u_noise , 2);

Ts_dec = t_noise_dec(2) - t_noise_dec(1);

% Dividiamo i dati di validazione e indentificazione

N = length(t_noise_dec);
N_30 = floor(N*0.3);

dati_val = iddata(y_noise_dec(1:N_30) , u_noise_dec(1:N_30) , Ts_dec);
dati_ident = iddata(y_noise_dec(N_30+1:end) , u_noise_dec(N_30+1:end) , Ts_dec);

% ARX con metodo AIC
ordini_ARX = [];
ord_max = 5;
for na = 1:ord_max
    for nb = 1:ord_max
        for ritardo = 1:ord_max
            ordini_ARX(na+nb+ritardo-2 , :) = [na , nb , ritardo]; 
        end
    end
end

V = arxstruc(dati_ident , dati_val , ordini_ARX);
ordine_modello_MDL = selstruc(V , 'MDL');
ordine_modello_AIC = selstruc(V , 'AIC');

%% Simulazione e Predizione
modello_arx_MDL = arx(dati_ident , ordine_modello_MDL);
modello_arx_AIC = arx(dati_ident , ordine_modello_AIC);
figure
compare(modello_arx_MDL , dati_val);
title("Simulazione MDL")
figure
compare(modello_arx_AIC , dati_val);
title("Simulazione AIC")

%% Analisi dei residui
figure
resid(dati_val , modello_arx_MDL);
title("Analisi Residui MDL")

figure
resid(dati_val , modello_arx_AIC);
title("Analisi Residui AIC")

% Il modello con MDL è quello migliore

%% Poli e Zeri
h = iopzplot(modello_arx_MDL);
showConfidence(h,3);

% La posizione tra i poli e zeri sono molto vicine, ma non dentro
% l'intervallo di confidenza

%% Modello OE
OE_cancellazione = oe(dati_ident , [3 1 1]); 
OE_senza_cancellazione = oe(dati_ident , [5 3 1]);

figure
resid(dati_val , OE_cancellazione);
title("Con Cancellazione");

figure
resid(dati_val , OE_senza_cancellazione);
title("Senza Cancellazione");

present(OE_senza_cancellazione)

figure
bodeplot(OE_senza_cancellazione)
hold on
bodeplot(FDT_non_parametrica)
%% Results

% Ghat; % model of exogenous term G(z)
% 
% Hhat; % model of noise term H(z)
% 
% 
