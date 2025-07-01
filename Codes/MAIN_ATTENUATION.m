clear
close all
clc

% dati

beta2=0.7;
xi=0.008681453698009;    % structural (short circuit)
wi=2*pi*27;            % short circuit
wi_hat=2*pi*27.711111;        % open circuit

Hsc=0.723547;
Hsh=0.122617;


% calcolo parametri

ki=sqrt((wi_hat^2-wi^2)/(wi^2));
ki_tilde=ki/sqrt(1-beta2);
K=ki*sqrt((beta2)/(1-beta2));

% calcolo attenuazione con formula

Adb_formula=20*log10( (ki_tilde^2+2*sqrt(2)*xi*sqrt(2+ki_tilde^2-2*K^2))/(4*xi*sqrt(1-xi^2)) );

% calcolo attenuazione con FRFs

Adb_experiments=20*log10(Hsc/Hsh);