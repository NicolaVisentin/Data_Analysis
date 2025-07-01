%% Cp1 non fornito

clear
close all
clc

% inserimento prime N frequenze proprie

omega_oc=[174.4631];    % open-circuit
omega_sc=[169.646];    % short-circuit

% inserimento coefficienti di accoppiamento k_i

k = sqrt((omega_oc.^2 - omega_sc.^2)./omega_sc.^2);

% calcolo del vettore di Cp

Cp = zeros(length(k),1);    % inizializzazione vettore di Cp
Cp(end) = 30.5;             % inserire il C_infinito [nF]
Cp(end)=Cp(end)*10^(-9);    % conversione Cp in [F]

for i = length(Cp)-1:-1:1
    Cp(i) = Cp(i+1) * (1 + k(i+1)^2);
end

% stampare il valore di Cp1 (prima mode)

fprintf(['Il Cp1 è ',num2str(Cp(1)*10^9),' nF\n'])

%% Cp1 fornito

clear
close all
clc

% define test parameters

Beta2=0.7; % [-]
Cp1=37.5;   % [nF]
f_oc=4.5927;  % open-circuit first natural frequency [Hz]
f_sc=4.3638;  % close-circuit first natural frequency [Hz]

% compute natural pulsations [rad/s] and convert capacity

w_oc=174.4631;
w_sc=169.646;

Cp1=Cp1*1e-9;
C2 = Cp1/Beta2;
% compute the equivalent capacity Ceq

Ceq=Cp1/(1-Beta2);

% compute the omega F parameter and tau_e

w_F=sqrt((w_sc.^2 + w_oc.^2)./2);
tau_e=1/w_F;

% compute the optimal Shunt resistance

R_sh=tau_e/Ceq;
fprintf(['The optimal Shunt resistance is ',num2str(R_sh),' [Ohm]\n'])
