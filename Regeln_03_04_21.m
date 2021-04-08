clear all; close all; clc;

%% Daten laden
load('Gruppe_C_Aufgabe.mat')
s = tf('s');

Mod = Modellfamilie;
% Zugriff auf ein Modell: Mod(:,:,1,1)
%% Normieren
e_max = 50;
u_max = 10;
d_max = 2;

% Erzeugung der Normierungsmatrizen
Du=u_max*eye(2);
De=e_max*eye(2);
Dd=d_max*eye(2);
%% Unsicherheitsbeschreibung

r0 = 10^(34.2/20);
rinf = 10^(35.2/20);
% r0 = 10^(48.3/20);
% rinf = 10^(49.3/20);
tau = 2;                % Ausprobieren
w_I = (tau * s + r0) / (tau / rinf * s + 1);
W_I = w_I * eye(2);

figure(3)
for i = 1:20
    sigma((inv(De)*Mod(:,:,i,1)*Du),'k')
    hold on
end
sigma(w_I,'magenta')
%% Parameter ranges angeben
% Den Zähler und den Nenner bekommen

num = zeros(length(Mod),4,3);           %Reihenfolge: 1-1,1-2,2-1,2-2
den = zeros(length(Mod),4,3);           %gleiche Reihenfolge
for i = 1:length(Mod)
    for j = 1:4
        [tempn,tempd] = tfdata(Mod(:,:,i,1));
        num(i,j,4-length(cell2mat(tempn(j))):end) = cell2mat(tempn(j));
        den(i,j,4-length(cell2mat(tempd(j))):end) = cell2mat(tempd(j));
    end
end

% Nominelle Werte bestimmen
num11 = reshape(mean(num(:,1,:)),[1,3]);
num12 = reshape(mean(num(:,2,:)),[1,3]);
num21 = reshape(mean(num(:,3,:)),[1,3]);
num22 = reshape(mean(num(:,4,:)),[1,3]);

den11 = reshape(mean(den(:,1,:)),[1,3]);
den12 = reshape(mean(den(:,2,:)),[1,3]);
den21 = reshape(mean(den(:,3,:)),[1,3]);
den22 = reshape(mean(den(:,4,:)),[1,3]);
% unsichere Parameter Ureal und nominelle Paramater festlegen
a1_11 = num11(1);
a2_11 = ureal('a2_11',num11(2),'Range',[min(num(:,1,2)) max(num(:,1,2))]);
a3_11 = ureal('a3_11',num11(3),'Range',[min(num(:,1,3)) max(num(:,1,3))]);

a1_21 = num21(1);
a2_21 = num21(2);
a3_21 = ureal('a3_21',num21(3),'Range',[min(num(:,3,3)) max(num(:,3,3))]);

a1_12 = num12(1);
a2_12 = ureal('a2_21',num12(2),'Range',[min(num(:,2,2)) max(num(:,2,2))]);
a3_12 = ureal('a3_12',num12(3),'Range',[min(num(:,2,3)) max(num(:,2,3))]);

a1_22 = num22(1);
a2_22 = num22(2);
a3_22 = ureal('a3_22',num22(3),'Range',[min(num(:,4,3)) max(num(:,4,3))]);

%Range für Nennerparameter
b1_11 = den11(1);
b2_11 = ureal('b2_11',den11(2),'Range',[min(den(:,1,2)) max(den(:,1,2))]);
b3_11 = ureal('b3_11',den11(3),'Range',[min(den(:,1,3)) max(den(:,1,3))]);

b1_21 = den21(1);
b2_21 = den21(2);
b3_21 = ureal('b3_21',den21(3),'Range',[min(den(:,3,3)) max(den(:,3,3))]);

b1_12 = den12(1);
b2_12 = den12(2);
b3_12 = ureal('b3_12',den12(3),'Range',[min(den(:,2,3)) max(den(:,2,3))]);
 
b1_22 = den22(1);
b2_22 = ureal('b2_22',den22(2),'Range',[min(den(:,4,2)) max(den(:,4,2))]);
b3_22 = ureal('b3_22',den22(3),'Range',[min(den(:,4,3)) max(den(:,4,3))]);
 
% Übertragungsfunktionen aufstellen
Gp11 = tf([a2_11 a3_11],[b1_11 b2_11 b3_11]);
Gp12 = tf([a2_12 a3_12],[b2_12 b3_12]);
Gp21 = tf([a2_21 a3_21],[b2_21 b3_21]);
Gp22 = tf([a2_22 a3_22],[b1_22 b2_22 b3_22]);

Gp_phys =uss([Gp11 Gp21; Gp12 Gp22]); % phys Modell
Gdp = ss([Gp11 Gp21; Gp12 Gp22]);% Störbertragungsfunktion
Gp=(inv(De)*Gp_phys*Du); % normiertes unsicheres Modell

Gn = tf(Gp.NominalValue); % normiertes nominelles Modell

%Sensitivität mit Performance-Gewicht W_P
A = 0.01; %Wir setzen das bei 0.01, da wir somit die komplementäre Sensitivität in eine stationäre Genauigkeit von +- 1% zwingen, da S+T = 1 (Gamma ist größer als eins, also wird die Anforderung +-3% erfüllt)
M = 2;
omega_B = 450;

W_P = eye(2)*(s/M+omega_B)/(s+omega_B*A);
W_U = eye(2);
Gd = eye(2);
%% Hinf
% Verallgemeinerte Strecke mit unstrukturierter Modellunsicherheit für
% Reglerauslegung
systemnames		= 'Gn Gd W_P W_I W_U';
inputvar			= '[uDel(2); r(2); d(2); u(2)]';
outputvar			= '[W_I; W_P; W_U; r-Gn-Gd]';
input_to_Gn		= '[u+uDel]';
input_to_Gd		= '[d]';
input_to_W_P		= '[r-Gn-Gd]';
input_to_W_I		= '[u]';
input_to_W_U		= '[u]';
sysoutname		= 'P';
cleanupsysic	= 'yes';
sysic;
P = minreal(P);

% Erstellung der unstrukturierten Modellunsicherheit
Delta = ultidyn('Delta',[2 2]);

% upper LFT
Pd = lft(Delta,P); %P Delta
Pd = minreal(Pd);

% Reglerauslegung
figure
[K1,CL1,GAM1,INFO1] = hinfsyn(Pd,2,2);
K1 = minreal(K1);
S1 = inv(eye(2) + Gn * K1); % Ausgangsseitige Sensitivität
T1 = eye(2) - S1;

% sigma(T1)
hold on
sigma(S1)
sigma(1/W_P)

legend('S1','1/W_P')
% legend('T1','S1','1/W_P')
N = lft(P, K1)
%% Dksyn
 %Sensitivität mit Performance-Gewicht W_P
A = 0.01; %Wir setzen das bei 0.01, da wir somit die komplementäre Sensitivität in eine stationäre Genauigkeit von +- 1% zwingen, da S+T = 1 (Gamma ist größer als eins, also wird die Anforderung +-3% erfüllt)
M = 2;
omega_B = 450;

W_P = eye(2)*(s/M+omega_B)/(s+omega_B*A);
W_U = eye(2);

Gd = eye(2);
% Verallgemeinerte Strecke mit unstrukturierter Modellunsicherheit für
% Reglerauslegung
systemnames		= 'Gn Gd W_P W_I W_U';
inputvar			= '[uDel(2); r(2); d(2); u(2)]';
outputvar			= '[W_I; W_P; W_U; r-Gn-Gd]';
input_to_Gn		= '[u+uDel]';
input_to_Gd		= '[d]';
input_to_W_P		= '[r-Gn-Gd]';
input_to_W_I		= '[u]';
input_to_W_U		= '[u]';
sysoutname		= 'P';
cleanupsysic	= 'yes';
sysic;
P = minreal(P);
% Erstellung der unstrukturierten Modellunsicherheit
Delta = ultidyn('Delta',[2 2]);

% upper LFT
Pd = lft(Delta,P); %P Delta
Pd = minreal(Pd);

% Reglerauslegung
figure
[K2,CL2,GAM2,INFO2] = dksyn(Pd,2,2);
K2 = minreal(K2);
S2 = inv(eye(2) + Gn * K2); % Ausgangsseitige Sensitivität
T2 = eye(2) - S2;

sigma(T2)
hold on
sigma(S2)
sigma(1/W_P)

legend('T2','S2','1/W_P')




%% 
% % Erstellung der unstrukturierten Modellunsicherheit
% Delta = ultidyn('Delta',[2 2]);
% 
% % LFT
% Pp=lft(Delta,P);
% Pp = minreal(Pp);
% 
% %Reglersynthese
% [K2,CL2,GAM2,INFO2] = dksyn(Pp,2,2);
% K2 = minreal(K2);
% S2 = inv(eye(2) + G * K2); % Ausgangsseitige Sensitivität
% T2 = eye(2) - S2;
% SI2 = inv(eye(2) + K2 * G); % Eingangsseitige Sensitivität
% TI2 = eye(2) - SI2;
%% Reglersynthese mit W_T
M = 1.03;
A = 0.01;
omega_T = 1000;

w_T = (s+omega_T/M)/(A*s+omega_T);
W_T = eye(2) * w_T;
W_U = eye(2);
Gd = eye(2);
%% DKsyn mit WT
% Verallgemeinerte Strecke mit unstrukturierter Modellunsicherheit für
% Reglerauslegung
systemnames		= 'Gn Gd W_T W_I W_U';
inputvar			= '[uDel(2); r(2); d(2); u(2)]';
outputvar			= '[W_I; W_T; W_U; r-Gn-Gd]';
input_to_Gn		= '[u+uDel]';
input_to_Gd		= '[d]';
input_to_W_T		= '[Gn+Gd]';
input_to_W_I		= '[u]';
input_to_W_U		= '[u]';
sysoutname		= 'P';
cleanupsysic	= 'yes';
sysic;
P = minreal(P);

% Erstellung der unstrukturierten Modellunsicherheit
Delta = ultidyn('Delta',[2 2]);

% upper LFT
Pd = lft(Delta,P); %P Delta
Pd = minreal(Pd);

% Reglerauslegung
figure
[K1,CL1,GAM1,INFO1] = dksyn(Pd,2,2);
K1 = minreal(K1);
S1 = inv(eye(2) + Gn * K1); % Ausgangsseitige Sensitivität
T1 = eye(2) - S1;

sigma(T1)
hold on
sigma(S1)
sigma(1/W_T)

legend('T1','S1','1/W_T')


% % Erstellung der unstrukturierten Modellunsicherheit
% Delta = ultidyn('Delta',[2 2]);
% 
% % LFT
% Pp=lft(Delta,P);
% Pp = minreal(Pp);
% 
% %Reglersynthese
% [K2,CL2,GAM2,INFO2] = dksyn(Pp,2,2);
% K2 = minreal(K2);
% S2 = inv(eye(2) + G * K2); % Ausgangsseitige Sensitivität
% T2 = eye(2) - S2;
% SI2 = inv(eye(2) + K2 * G); % Eingangsseitige Sensitivität
% TI2 = eye(2) - SI2;



