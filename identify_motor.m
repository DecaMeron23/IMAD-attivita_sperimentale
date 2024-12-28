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
% Valutazione del tempo di campionamento
Ts = t_raw_step(2) - t_raw_step(1);

% Valutazione del ritardo k, si osserva che il ritardo è di circa 5 passi
figure
hold on
plot(t_raw_step , u_raw_step ,"DisplayName","PWM [$\%$]");
plot(t_raw_step , y_raw_step , "DisplayName", "rpm $\left[ \frac{{giri}}{{min}}\right]$");
title("Valutazione Ritardo");
xlabel("Tempo [$s$]" , Interpreter="latex");
xlim([11.0640 11.1290])
ylim([0.0 54.9])
rpmleftfracgiriminrightLine = findobj(gcf, "DisplayName", "rpm $\left[ \frac{{giri}}{{min}}\right]$");
datatip(rpmleftfracgiriminrightLine,11.11,20.23);
PWMLine = findobj(gcf, "DisplayName", "PWM [$\%$]");
datatip(PWMLine,11.08,50);
legend(Location="southeast" , Interpreter="latex");

% Si osserva un ritardo puro di 5 passi
ritardo = 5;

% Plot dello scalino con rumore
figure
hold on
plot(t_raw_noise , u_raw_noise , "DisplayName" , "PWM [$\%$]");
plot(t_raw_noise, y_raw_noise , "DisplayName" , "rpm $\left[ \frac{{giri}}{{min}}\right]$");
title("Scalino con Rumore");
xlabel("Tempo [$s$]" , Interpreter="latex");
legend(Location="southeast");

%% Pulizia segnali

% Si è scelto l'intervallo tra 6s e 23s per avere un segnale stazionario
intervallo = floor(6/Ts) : floor(23/Ts);

t_noise = t_raw_noise(intervallo); 
u_noise = u_raw_noise(intervallo);
y_noise = y_raw_noise(intervallo);

% Rimuoviamo la media
u_media = mean(u_noise);
u_noise = u_noise - u_media;

y_media = mean(y_noise);
y_noise = y_noise - y_media;

% Riscaliamo il tempo in modo che parta da 0, dopo aver tolto i secondi
% iniziali (riferimento riga 66)
t_noise = t_noise - t_noise(1); 

% plot dei segnali elaborati
figure
hold on
plot(t_noise , u_noise , DisplayName="PWM [$\%$]");
plot(t_noise , y_noise , DisplayName="rpm $\left[ \frac{{giri}}{{min}}\right]$");
title("Segnali Elaborati");
xlabel("Tempo [$s$]" , Interpreter="latex");
legend();

%% Identificazione non parametrica

% Determinamo la finestra di lavoro
fs = 1/Ts;
ws = 0.05;
w = 0:ws:fs*pi;

data = iddata(y_noise, u_noise, Ts);

FDT_non_parametrica = spafdr(data , [] , w);

figure
h = bodeplot(FDT_non_parametrica);
showConfidence(h,3)

%% Identificazione parametrica
% Dal diagramma di bode del modulo osserviamo una pulsazione critica pari
% a circa 16 rad/s che corrisponde a circa 2.55hz quindi ricampioniamo il segnale
% a circa 50hz (10 w_c/2pi < fs < 30 w_c/2pi)

% Ricampioniamo a 50hz,corrispondente alla riduzione di un fattore 4 ->
% 50Hz = (1/0.005s)/4 = 200Hz/4
t_noise_dec = decimate(t_noise , 4);
y_noise_dec = decimate(y_noise , 4);
u_noise_dec = decimate(u_noise , 4);

% Estraiamo il nuovo tempo di campionamento
Ts_dec = t_noise_dec(2) - t_noise_dec(1);

% Dividiamo i dati in validazione e identificazione
N = length(t_noise_dec);
N_30 = floor(N*0.3);

dati_val = iddata(y_noise_dec(1:N_30) , u_noise_dec(1:N_30) , Ts_dec);
dati_ident = iddata(y_noise_dec(N_30+1:end) , u_noise_dec(N_30+1:end) , Ts_dec);

%% Ricerca del Modello ARX migliore tramite formule di complessità
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
ordine_modello = selstruc(V , 'AIC');
% Tramite la formula di complessità si è notato che il ritardo puro che dà
% il modello migliore è pari a 1 e non a 5 (come osservato in precedenza).
% Si sono svolte delle prove con un ritardo puro a 5 ma il modello migliore
% trovato tramite la formula di complessità risultava essere sempre
% peggiore del modello utilizzato da qui in avanti (ovvero quello con
% ritardo puro pari a 1)

%% Analisi del modello
modello_arx = arx(dati_ident , ordine_modello);

% Simulazione
figure
compare(modello_arx , dati_val);

% Analisi dei residui
figure
resid(dati_val , modello_arx);

%% Verifica Cancellazi tra poli e zeri
h = iopzplot(modello_arx);
showConfidence(h,3);

% Non c'è nessuna possibile cancellazione

%% Modello OE
OE = oe(dati_ident , ordine_modello);

% Analisi dei Residui
figure
resid(dati_val , OE);
% Si osserva che la crosscorrelazione è molto buona, mentre
% autocorrelazione non lo è (come ci aspettavamo).

% Verifica se si possono effettuare delle riduzioni
present(OE)

% Comparazione del diagramma di bode con la stima non parametrica, si
% osserva la parte esogena viene modellata bene.
figure
bodeplot(OE);
hold on
bodeplot(FDT_non_parametrica);
xlim([1e-1 1e2]);
legend(["OE" , "FdT non parametrica"]);

% Analisi dello spettro, come ci si aspettava OE non modella questa parte
figure
spectrumplot(OE);
hold on
spectrumplot(FDT_non_parametrica);
xlim([1e-1 1e2]);
legend(["OE" , "FdT non parametrica"]);


%% Modellazione del rumore

% Stimiamo il rumore
y_sim_OE = sim(OE, dati_val);
errore_OE = dati_val.OutputData - y_sim_OE.OutputData;
v_spec = spafdr(iddata(errore_OE));

% plot dello spettro
figure
spectrumplot(v_spec);

% Dall'analisi dello spettro si ipotizza di utilizzare 2 poli e 2 zeri
ordini_BJ = [ordine_modello(2) , 2, 2, ordine_modello(1) , ordine_modello(3)];

% Modello Box-Jenkins
BJ = bj(dati_ident , ordini_BJ);

% Verifica se si possono effettuare delle riduzioni
present(BJ)

% Analisi dei residui
figure
resid(dati_val , BJ);

% In generale sia la parte esogena che la parte del disturbo vengono
% modellate bene, seppur non ottimali.

%% Confronto con stima non parametrica

% Plot dei diagrammi di bode per il confronto con la stima non parametrica
figure
bodeplot(BJ);
hold on
bodeplot(FDT_non_parametrica)

legend(["BJ" , "FdT non parametrica"]);

% Plot dello spettro del rumore per il confronto con la stima non 
% parametrica
figure
spectrumplot(BJ);
hold on
spectrumplot(FDT_non_parametrica)

legend(["BJ" , "FdT non parametrica"]);

%% Simulazione dei tre modelli

% Simulazione ARX
y_arx = sim(modello_arx, dati_val);
figure
compare(modello_arx , dati_val);

% Simulazione OE
y_OE = sim(OE , dati_val);
figure
compare(OE, dati_val);

% Simulazione BJ
y_BJ = sim(BJ , dati_val);
figure
compare(BJ , dati_val);

% Dalle simulazioni il Box-Jenkins è il migliore, lo riaddestriamo
% sull'intero data set

%% BJ data set completo

dati = iddata(y_noise_dec , u_noise_dec , Ts_dec);

BJ = bj(dati , ordini_BJ);

BJ

%% Results

% Definiamo Ghat
numerator_G = tf(BJ.B , 1 , Ts_dec);
numerator_G.Variable = "z^-1";
numerator_G.Denominator = 1;

denominator_G = tf(BJ.F , 1 , Ts_dec);
denominator_G.Variable = "z^-1";
denominator_G.Denominator = [0 1];

Ghat = numerator_G/denominator_G % model of exogenous term G(z)

% Definiamo Hhat
numerator_H = tf(BJ.C , 1 , Ts_dec);
numerator_H.Variable = "z^-1";
numerator_H.Denominator = 1;

denominator_H = tf(BJ.D , 1 , Ts_dec);
denominator_H.Variable = "z^-1";
denominator_H.Denominator = 1;

Hhat = numerator_H / denominator_H % model of noise term H(z)


%% Buon Anno e Buone Feste
% Dati per il triangolo (chioma dell'albero)
x1 = [-2 0 2]; % Base del triangolo
y1 = [0 3 0];  % Altezza del triangolo

x2 = [-1.5 0 1.5];
y2 = [1 4 1];

x3 = [-1 0 1];
y3 = [2 5 2];

% Dati per il tronco
x_tronco = [-0.5 0.5 0.5 -0.5];
y_tronco = [0 -0.5 -1 -1];

% Disegna l'albero
figure;
hold on;
fill(x1, y1, [0 0.5 0], 'EdgeColor', 'none'); % Primo triangolo (verde)
fill(x2, y2, [0 0.6 0], 'EdgeColor', 'none'); % Secondo triangolo (verde scuro)
fill(x3, y3, [0 0.7 0], 'EdgeColor', 'none'); % Terzo triangolo (verde scuro)
fill(x_tronco, y_tronco, [0.5 0.25 0], 'EdgeColor', 'none'); % Tronco (marrone)

% Aggiungi decorazioni
scatter(0, 5.2, 100, 'yellow', 'filled'); % Stella sulla punta
scatter(rand(1, 10)*4-2, rand(1, 10)*3, 50, 'red', 'filled'); % Palline rosse

% Migliora l'aspetto
axis equal;
axis([-3 3 -2 6]);
title('Buone Feste!!');
hold off;
