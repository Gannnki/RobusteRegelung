close all
clear all
clc

load('Gruppe_C_Aufgabe.mat')
s = tf('s');

% Normierung
umax = 10;
emax = 50 ;
dmax = 2;

%Größe der Modellfamilie
[INnum,OUTnum,MODnum] = size(Modellfamilie);

Du=umax*eye(2);
De=emax*eye(2);
Dd=dmax*eye(2);
fprintf('**************************\n')
fprintf('System Overview : ')
fprintf('The model family contains %d models.\n',MODnum)
fprintf('Each model has %d inputs.\n',INnum)
fprintf('               %d outputs.\n',OUTnum)
%% System characteristics : Aufgabe 1

for i = 1:MODnum
    G(:,:,i,1)  = Modellfamilie(:,:,i,1);
end
fprintf('**************************\n')
fprintf('System characteristics:\n')
for i = 1:MODnum
    G11 = G(1,1,i,1);% Input 1 to output 1
    G12 = G(1,2,i,1);% Input 1 to output 2
    G21 = G(2,1,i,1);% Input 2 to output 1
    G22 = G(2,2,i,1);% Input 2 to output 2
    if all(real(pole(G11)))&&all(real(pole(G22)))
        fprintf('The %d th. MIMO modell is asymptotically stable.\n',i)
    end
end
fprintf('**************************\n')
%% System uncertainty

% To diffrentiate from Assignment 1. we might need to call the for-loop
% again
% 4 Variables for G11 : G11_a G11_b G11_c G11_d
% 3 Variables for G12 : G12_a G12_b G12_c 
% 2 Variables for G21 : G21_a G21_b
% 4 Variables for G22 : G22_a G22_b G22_c G22_d

for i = 1:MODnum
    G11 = G(1,1,i,1);% Input 1 to output 1
    G12 = G(1,2,i,1);% Input 1 to output 2
    G21 = G(2,1,i,1);% Input 2 to output 1
    G22 = G(2,2,i,1);% Input 2 to output 2
    % G11
    G11_a(i) = get(G11).K;
    G11_b(i) = -cell2mat(get(G11).Z);
    P = -cell2mat(get(G11).P);
    G11_c(i) = P(1,:);
    G11_d(i) = P(2,:);
    % G12
    G12_a(i) = get(G12).K;%Gain
    G12_b(i) = -cell2mat(get(G12).P);%Zero
%     G12_c(i) = cell2mat(get(G12).P);%Pole
    % G21
    G21_a(i) = get(G21).K;
    G21_b(i) = -cell2mat(get(G21).Z);%Pole
    G21_c(i) = -cell2mat(get(G21).P);
    
    % G22
    G22_a(i) = get(G22).K;
    G22_b(i) = -cell2mat(get(G22).Z);
    P = -cell2mat(get(G22).P);
    G22_c(i) = P(1,:);
    G22_d(i) = P(2,:);
end
%Calculate mean value, maximal and minimal values
%G11
G11_a_m = mean(G11_a);G11_a_t = max(G11_a);G11_a_d = min(G11_a);% top and down
G11_b_m = mean(G11_b);G11_b_t = max(G11_b);G11_b_d = min(G11_b);
G11_c_m = mean(G11_c);G11_c_t = max(G11_c);G11_c_d = min(G11_c);
G11_d_m = mean(G11_d);G11_d_t = max(G11_d);G11_d_d = min(G11_d);
%G12
G12_a_m = mean(G12_a);G12_a_t = max(G12_a);G12_a_d = min(G12_a);
G12_b_m = mean(G12_b);G12_b_t = max(G12_b);G12_b_d = min(G12_b);
% G12_c_m = mean(G12_c);G12_c_t = max(G12_c);G12_c_d = min(G12_c);
%G21
G21_a_m = mean(G21_a);G21_a_t = max(G21_a);G21_a_d = min(G21_a);
G21_b_m = mean(G21_b);G21_b_t = max(G21_b);G21_b_d = min(G21_b);
G21_c_m = mean(G21_c);G21_c_t = max(G21_c);G21_c_d = min(G21_c);
%G22
G22_a_m = mean(G22_a);G22_a_t = max(G22_a);G22_a_d = min(G22_a);% top and down
G22_b_m = mean(G22_b);G22_b_t = max(G22_b);G22_b_d = min(G22_b);
G22_c_m = mean(G22_c);G22_c_t = max(G22_c);G22_c_d = min(G22_c);
G22_d_m = mean(G22_d);G22_d_t = max(G22_d);G22_d_d = min(G22_d);

% pertubated transfer function
%G11
G11_a_p = ureal('G11a', G11_a_m, 'Range', [G11_a_d G11_a_t]);
G11_b_p = ureal('G11b', G11_b_m, 'Range', [G11_b_d G11_b_t]);
G11_c_p = ureal('G11c', G11_c_m, 'Range', [G11_c_d G11_c_t]);
G11_d_p = ureal('G11d', G11_d_m, 'Range', [G11_d_d G11_d_t]);
%G12
G12_a_p = ureal('G12a', G12_a_m, 'Range', [G12_a_d G12_a_t]);
G12_b_p = ureal('G12b', G12_b_m, 'Range', [G12_b_d G12_b_t]);
% G12_c_p = ureal('G12c', G12_c_m, 'Range', [G12_c_d G12_c_t]);
%G21
G21_a_p = ureal('G21a', G21_a_m, 'Range', [G21_a_d G21_a_t]);
G21_b_p = ureal('G21b', G21_b_m, 'Range', [G21_b_d G21_b_t]);
G21_c_p = ureal('G21c', G21_c_m, 'Range', [G21_c_d G21_c_t]);
%G22
G22_a_p = ureal('G22a', G22_a_m, 'Range', [G22_a_d G22_a_t]);
G22_b_p = ureal('G22b', G22_b_m, 'Range', [G22_b_d G22_b_t]);
G22_c_p = ureal('G22c', G22_c_m, 'Range', [G22_c_d G22_c_t]);
G22_d_p = ureal('G22d', G22_d_m, 'Range', [G22_d_d G22_d_t]);

G11_p = G11_a_p*(s+G11_b_p)/((s+G11_c_p)*(s+G11_d_p));
G21_p = G21_a_p*(s+G21_b_p)/(s+G21_c_p);
G12_p = G12_a_p/(s+G12_b_p);
G22_p = G22_a_p*(s+G22_b_p)/((s+G22_c_p)*(s+G22_d_p));
G_norm = inv(De)*G*Du;
Gp = uss([G11_p,G12_p;G21_p,G22_p]);
Gp_norm = inv(De)*Gp*Du;
% k is the factor od disturbance
k=1;% We can choose this factor based on some literature recommendation
Gd = k*eye(2);
% What interesting is : in state space modell under perspective of
% practical usage the disturbance matrix is ofen chosen as diag(1e-1).(10^-1)
% But yet not clear whether it is directly applicable.
Gd = inv(De)*Gd*Dd;
% Also the second method is very sensitive to Gd --> reason? unknown. 
Gn = tf(Gp_norm.NominalValue);%Gpn
[N,Delta]=lftdata(Gp);%Verallgemeinerte Strecke
%% Gewicht > obere Schranke

% gleich wie Übung 12, nehmen wir eine 
% multiplikative, eingangsseitige Modellunsicherheit für die Strecke an.
% Die Formel steht auf Seite 69.(SISO) und 153.(MIMO)
% See Skript Page 151

figure(1)

for i =1:MODnum
    sigma(inv(Gn)*(G_norm(:,:,i,1)-Gn),'r-.');% Abweichung
    drawnow; 
    hold on
    grid on
end

r0 = db2mag(-10);%See script from Prof.King Page 84.
rinf = db2mag(-16);
tau = 0.004;
wI = (tau*s+r0)/(tau/rinf*s+1);
sigma(wI,'g')
xlim([1e-4 1e5])
hold off
%% Requirements on controller

% 1 Führungsgröße : schnell w_b > 0.025 rad/s --> untere Grenze der Bandbreite
% 2 bleibende Regelfehler : e < 3 % > 1% --> A = 0.01
% 3 Änderungsgeschwindigkeit der Stellgröße U < 1256.64 A/sec
% Page 39. eine Foderung für KS mit Gewicht von 'Wu'
% Störung : disturbance rejection > Anforderung an S > large loop gain
%% Definition der Gewichte

Wu = tf(1)*eye(2);
Wb = 0.025;% 40
M = 1.5;
A = 0.01 ;%wie in Übung 4
Wp = tf([1/M Wb],[1 A*Wb])*eye(2);
Wi = wI*eye(2);
%% H-inf Synthese

systemnames		= 'Gn Gd Wp Wi Wu';
inputvar			= '[uDel(2); r(2); d(2); u(2)]';
outputvar			= '[Wi; Wp; Wu; r-Gn-Gd]';
input_to_Gn		= '[u+uDel]';
input_to_Gd		= '[d]';
input_to_Wp		= '[r-Gn-Gd]';
input_to_Wi		= '[u]';
input_to_Wu		= '[u]';
sysoutname		= 'P';
cleanupsysic	= 'yes';
sysic;

P1 = minreal(P);

% Erstellung der unstrukturierten Modellunsicherheit
Delta = ultidyn('Delta',[2 2]);

%LFT
Pp=lft(Delta,P1);
Pp = minreal(Pp);

%Reglersynthese
[K1,CL1,GAM1,INFO1] = hinfsyn(Pp,2,2);
K1 = minreal(K1);
S1 = inv(eye(2) + Gn * K1);
T1 = eye(2) - S1;
SI1 = inv(eye(2) + K1 * Gn);
TI1 = eye(2) - SI1;

% Auswertung (tbc)
figure(200)
sigma(S1)
hold on
sigma(inv(Wp))
grid on
xlim([1e-4 1e5])
title('Sensitivity')

figure(201)
sigma(K1*S1)
hold on
sigma(inv(Wu))
hold off
grid on
xlim([1e-4 1e5])
title('KS')

figure(202)
sigma(TI1)
hold on
sigma(inv(wI))
grid on
xlim([1e-4 1e5])
title('TI')
%% Mu-Synthese

[K2,CL2] = musyn(Pp,2,2);% gleiche Systemstruktur
K2  = minreal(K2);
S2  = inv(eye(2) + Gn * K2); % Ausgangsseitige Sensitivität
T2  = eye(2) - S2;
SI2 = inv(eye(2) + K2 * Gn); % Eingangsseitige Sensitivität
TI2 = eye(2) - SI2;% Gamma = 0.945 = Gam1
% Gam is normally smaller than 1, if Gamma bigger than 3,
% it is ofen assumed that the adjustment of W1, W2 and W3 is necessary.
%% Comparison and Evaluation Hinf

% nominelle Stabilität
% See Übung 8 Lösung Code line 11
if all(pole(T1))
    fprintf('**************************\n')
    fprintf('The Controller using H_inf Synthesis has Nominal Stability.\n')
    fprintf('**************************\n')
else
    fprintf('**************************\n')
    fprintf('The Controller using H_inf Synthesis does not have Nominal Stability.\n')
    fprintf('**************************\n')
end
if all(pole(T2))
    fprintf('**************************\n')
    fprintf('The Controller using Miu Synthesis with DK-Iteration has Nominal Stability.\n')
    fprintf('**************************\n')
else
    fprintf('**************************\n')
    fprintf('The Controller using Miu Synthesis with DK-Iteration does not have Nominal Stability.\n')
    fprintf('**************************\n')
end

[Plant,Delta,BLKRS] = lftdata(Pp); % Verallg. Strecke, Delta, Block für RS

% N-Delta-Struktur / Regler K1 einsetzen
N1 = lft(Plant,K1);        % lower LFT
N1.OutputName;
N1.InputName;

% Frequency Response Data für mussv
omega   = logspace(-5,4,100);% 1000 also tried, sadly too slow
Nfrd 	= frd(N1,omega); % RP: N-Delta-Struktur
Mfrd    = Nfrd(1:2,1:2); % RS: M-Delta-Struktur (N11)
N22frd  = Nfrd(3:6,3:6); % NP: Nur WP*S / WU*u (5:6) wird nicht betrachtet

% Frequency Response Data für mussv
omega   = logspace(-4,2,100);
Nfrd 	= frd(N1,omega); % RP: N-Delta-Struktur
Mfrd    = Nfrd(1:2,1:2); % RS: M-Delta-Struktur (N11)
N22frd  = Nfrd(3:6,3:6); % NP: Nur WP*S / WU*u (5:6) wird nicht betrachtet

% Unsicherheitsblöcke
BLKNP = [4 4];            % vollbesetzte 4x4-Matrix für Untersuchung der NP
BLKRP = [BLKRS.Size;4 4]; % Block für RP (2 für WP, 2 für Wu);  Erst BLKRS

% mussv
boundsNP = mussv(N22frd,BLKNP);
boundsRS = mussv(Mfrd,BLKRS);
boundsRP = mussv(Nfrd,BLKRP);

figure(300)
subplot(3,1,1)
semilogx(boundsNP,'k')
ylabel('NP Hinf')
hold on
grid on
xlim([1e-4 1e2])
subplot(3,1,2)
semilogx(boundsRS,'k')
ylabel('RS Hinf')
hold on
grid on
xlim([1e-4 1e2])
subplot(3,1,3)
semilogx(boundsRP,'k')
ylabel('RP Hinf')
hold on
grid on
xlim([1e-4 1e2])

% -- robustperf und robuststab
F1 = lft(Pp,K1);
[PERFMARG1,PERFMARGUNC,REPORT,INFO] = robustperf(F1);
subplot(3,1,3)
semilogx(INFO.MussvBnds(:,1),'r--')
[STABMARG1,DESTABUNC,REPORT,INFO] = robuststab(F1);
subplot(3,1,2)
semilogx(INFO.MussvBnds(:,1),'r--')
%% Comparison and Evaluation Miu

% N-Delta-Struktur / Regler K1 einsetzen
N2 = lft(Plant,K2);        % lower LFT
N2.OutputName;
N2.InputName;

% Frequency Response Data für mussv
omega   = logspace(-5,4,100);% 1000 also tried, sadly too slow
Nfrd 	= frd(N2,omega); % RP: N-Delta-Struktur
Mfrd    = Nfrd(1:2,1:2); % RS: M-Delta-Struktur (N11)
N22frd  = Nfrd(3:6,3:6); % NP: Nur WP*S / WU*u (5:6) wird nicht betrachtet

% Formel
% F = N22+N21*Delta*(I-N11*Delta)^(-1)*N12;
% Übertragungsmatrix des CLs

% Unsicherheitsblöcke

BLKNP = [4 4];            % vollbesetzte 4x4-Matrix für Untersuchung der NP
BLKRP = [BLKRS.Size;4 4]; % Block für RP (2 für WP, 2 für Wu);  Erst BLKRS

% mussv
boundsNP = mussv(N22frd,BLKNP);
boundsRS = mussv(Mfrd,BLKRS);
boundsRP = mussv(Nfrd,BLKRP);

figure(301)
subplot(3,1,1)
semilogx(boundsNP,'k')
ylabel('NP Miu')
grid on
hold on
xlim([1e-4 1e5])
subplot(3,1,2)
semilogx(boundsRS,'k')
ylabel('RS Miu')
grid on
hold on
xlim([1e-4 1e5])
subplot(3,1,3)
semilogx(boundsRP,'k')
ylabel('RP Miu')
grid on
hold on
xlim([1e-4 1e5])

% -- robustperf und robuststab
F2 = lft(Pp,K2);
[PERFMARG2,PERFMARGUNC,REPORT,INFO] = robustperf(F2);
subplot(3,1,3)
semilogx(INFO.MussvBnds(:,1),'r--')
[STABMARG2,DESTABUNC,REPORT,INFO] = robuststab(F2);
subplot(3,1,2)
semilogx(INFO.MussvBnds(:,1),'r--')
hold on
%% Modellreduktion

% K1 Ordnung 8; K2 Ordnung 22.
% Wählen wir den zweite Entwurf zur Modellreduktion an. ;)
% Die Vorgehensweise ist ähnlich wie Übung 
n_des = 10;
% a b c 3 Variante
K2a = reduce(K2,n_des);
S2a = inv(eye(2) + Gn * K2a);
T2a = eye(2) - S2a;

K2b = balmr(K2,1,n_des);% Square-root balanced truncation
S2b = inv(eye(2) + Gn * K2b);
T2b = eye(2) - S2b;

K2c = ncfmr(K2,n_des); % balanced truncation model reduction for normalized G
S2c = inv(eye(2) + Gn * K2c);
T2c = eye(2) - S2c;

figure(400)
sigma(K2,K2a,K2b,K2c)
grid on
xlim([1e-4 1e5])
legend('origin','reduce','balmr','ncfmr')
%% Comparison and Evaluation Miu with model reduction

% N-Delta-Struktur / Regler K1 einsetzen
N2 = lft(Plant,K2);        % lower LFT
N2_a =  lft(Plant,K2a);
N2_b =  lft(Plant,K2b);
N2_c =  lft(Plant,K2c);


% Frequency Response Data für mussv
omega   = logspace(-5,4,100);% 1000 also tried, sadly too slow
Nfrd 	= frd(N2,omega); % RP: N-Delta-Struktur
Nfrd_a 	= frd(N2_a,omega);
Nfrd_b 	= frd(N2_b,omega);
Nfrd_c 	= frd(N2_c,omega);

Mfrd    = Nfrd(1:2,1:2); % RS: M-Delta-Struktur (N11)
Mfrd_a  = Nfrd_a(1:2,1:2);
Mfrd_b  = Nfrd_b(1:2,1:2);
Mfrd_c  = Nfrd_c(1:2,1:2);
N22frd  = Nfrd(3:6,3:6); % NP: Nur WP*S / WU*u (5:6) wird nicht betrachtet
N22frd_a  = Nfrd_a(3:6,3:6);
N22frd_b  = Nfrd_b(3:6,3:6);
N22frd_c  = Nfrd_c(3:6,3:6);

% Unsicherheitsblöcke

BLKNP = [4 4];            % vollbesetzte 4x4-Matrix für Untersuchung der NP
BLKRP = [BLKRS.Size;4 4]; % Block für RP (2 für WP, 2 für Wu);  Erst BLKRS

% mussv
boundsNP = mussv(N22frd,BLKNP);
boundsNP_a = mussv(N22frd_a,BLKNP);
boundsNP_b = mussv(N22frd_b,BLKNP);
boundsNP_c = mussv(N22frd_c,BLKNP);

boundsRS = mussv(Mfrd,BLKRS);
boundsRS_a = mussv(Mfrd_a,BLKRS);
boundsRS_b = mussv(Mfrd_b,BLKRS);
boundsRS_c = mussv(Mfrd_c,BLKRS);

boundsRP = mussv(Nfrd,BLKRP);
boundsRP_a = mussv(Nfrd_a,BLKRP);
boundsRP_b = mussv(Nfrd_b,BLKRP);
boundsRP_c = mussv(Nfrd_c,BLKRP);

figure(401)
subplot(3,1,1)
semilogx(boundsNP,'y')
ylabel('NP Miu')
grid on
hold on
semilogx(boundsNP_a,'m')
semilogx(boundsNP_b,'c')
semilogx(boundsNP_c,'r')
legend('origin','reduce','balmr','ncfmr')
xlim([1e-4 1e5])
hold off


subplot(3,1,2)
semilogx(boundsRS,'y')
ylabel('RS Miu')
grid on
hold on
semilogx(boundsRS_a,'m')
semilogx(boundsRS_b,'c')
semilogx(boundsRS_c,'r')
legend('origin','reduce','balmr','ncfmr')
xlim([1e-4 1e5])
hold off

subplot(3,1,3)
semilogx(boundsRP,'y')
ylabel('RP Miu')
grid on
hold on
semilogx(boundsRP_a,'m')
semilogx(boundsRP_b,'c')
semilogx(boundsRP_c,'r')
legend('origin','reduce','balmr','ncfmr')
xlim([1e-4 1e5])
hold off
suptitle('Comparison between original controller and reduced')
% Die Kurven liegen sehr nah aneinander, deswegen spielt die
% Modellreduktion fast keine Rolle bei Systemeigenschaften.

%% Führungs -und Störgrößenverhalten

% Führungsgrößensprünge
r1 = [-1; 1];
r2 = [1; -1];% Unterschiedliche Anregungsmoden

%  Störgrößensprünge
d1 = [-1; 1];
d2 = [1; -1];% Unterschiedliche Anregungsmoden

figure(500)
subplot(1,2,1)
step(T1*r1,'k',T2*r1,'b')
grid on
title('y')
subplot(1,2,2)
step(K1*S1*r1,'k',K2*S2*r1,'b')
title('u')
grid on
suptitle('Führungsgrößenverhalten unter r1 = [-1;1] norminell')

figure(501)
subplot(1,2,1)
step(S1*Gd*d1,'k',S2*Gd*d1,'b')
grid on
title('y')
subplot(1,2,2)
step(-K1*S1*Gd*d1,'k',-K2*S2*Gd*d1,'b')
title('u')
grid on
suptitle('Störgrößenverhalten unter d1 = [-1;1] norminell')


figure(502)
subplot(1,2,1)
step(T1*r2,'k',T2*r2,'b')
grid on
title('y')
subplot(1,2,2)
step(K1*S1*r2,'k',K2*S2*r2,'b')
title('u')
grid on
suptitle('Führungsgrößenverhalten unter r2 = [1;-1] norminell')

figure(503)
subplot(1,2,1)
step(S1*Gd*d2,'k',S2*Gd*d2,'b')
grid on
title('y')
subplot(1,2,2)
step(-K1*S1*Gd*d2,'k',-K2*S2*Gd*d2,'b')
title('u')
grid on
suptitle('Störgrößenverhalten unter d2 = [1;-1] norminell')
%% Was noch dazu (at WEEKEND ! > bevor 12.April )

% 1 : eine genaue Beschreibung des Gds > zurzeit nur als Einheitsmatrix
% angenommen

% 2 : nominelle Stabilität (SISO Krieterien noch gültig?)

% 3 : eine richtige Intepretation zur Aufgabe 5 

% 4 : der Befehl dysyn oder musyn; obsolete warning by Dysyn and musyn
% generates an overly high order Controller which made a model reduction 

% 5 : change plot number 

% 6 : color is completely wrong in figure(302)