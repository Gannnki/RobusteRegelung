%% Daten laden
load('Gruppe_C_Aufgabe.mat')
s = tf('s');

Mod = Modellfamilie;

% Zugriff auf ein Modell: Mod(:,:,1,1)

%% Systemeigenschaften --> bisher kenne ich nur Stabilität? zählt Performance auch dazu?

% Stabilität prüfen --> Nyquist
% nyquist(Mod)
% alle sehen ziemlich gleich aus

% # instabile Pole
m_0 = 0;
% # Integratorpole
l_0 = 0;

delta_phi = m_0*pi + l_0*pi/2;
% Keins der Modelle hat instabile Pole oder Integratorpole --> delta_phi
% sollte immer = 0 sein

% Dann müsste 2 zu 2 instabil sein?

%% Nominelles Modell aufstellen
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

% Die könnte ich dann alle addieren und dann mitteln und maximale
% Abweichung berechnen --> nominelles Modell
num11 = reshape(mean(num(:,1,:)),[1,3]);
num12 = reshape(mean(num(:,2,:)),[1,3]);
num21 = reshape(mean(num(:,3,:)),[1,3]);
num22 = reshape(mean(num(:,4,:)),[1,3]);

den11 = reshape(mean(den(:,1,:)),[1,3]);
den12 = reshape(mean(den(:,2,:)),[1,3]);
den21 = reshape(mean(den(:,3,:)),[1,3]);
den22 = reshape(mean(den(:,4,:)),[1,3]);

G_n11 = tf(num11,den11);
G_n12 = tf(num12,den12);
G_n21 = tf(num21,den21);
G_n22 = tf(num22,den22);

G_n = [G_n11 G_n21; G_n12 G_n22];

%% Normieren
e_max = 50;
u_max = 10;
d_max = 2;

% Erzeugung der Normierungsmatrizen
Du=u_max*eye(2);
De=e_max*eye(2);
Dd=d_max*eye(2);

% Normierung der Übertragungsfunktionen
Gn_norm=(inv(De)*G_n*Du);


figure(2)
for i = 1:20
    sigma((inv(De)*Mod(:,:,i,1)*Du),'k')
    hold on
end
sigma(Gn_norm,'r')

%% Unsicherheitsbeschreibung

r0 = 10^(34.2/20);
rinf = 10^(35.2/20);
% r0 = 10^(48.3/20);
% rinf = 10^(49.3/20);
tau = 2;                % Ausprobieren
w_m = (tau * s + r0) / (tau / rinf * s + 1);

figure(3)
for i = 1:20
    sigma((inv(De)*Mod(:,:,i,1)*Du),'k')
    hold on
end
sigma(w_m,'magenta')

G_p = Gn_norm.*(eye(2) + w_m.*ones(2,2));




% umwandeln in ein Modell, das weiter genutzt werden kann --> was sind die
% Kriterien?


% [u1_phys; u2_phys] = [10 0;10 0] .*[u1; u2];


