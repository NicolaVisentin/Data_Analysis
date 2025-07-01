clear
close all
clc

%% PLOTS CONTROL

plot_time=0;        % signals in time domain
plot_freq=0;        % signals in frequency domain
plot_powcross=1;    % power and cross spectra
plot_coh=0;         % coherence function
plot_tf=0;          % transfer function estimation (H1 and H2)
plot_modrec=0;      % reconstruction of estimated transf function

%% TRANSFER FUNCTION ESTIMATION FROM EXPERIMENTAL DATA

%% Load experimental data

% load data: we expect 2 records: random input (white noise) (first in the 
% data file) and its output (second in the data file); we also need to know
% the sampling frequency and the sensitivity.

[data_file,path]=uigetfile;
cd(path);
load(data_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO REMOVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "sistema" il file "di prova" per renderlo simile a quello che avremo

%Data=Data(:,1:2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set/compute test parameters

% change names to the data so that the code runs

x=Dati(:,2);    % input is called x
y=Dati(:,1);    % output is called y
fsamp=2048;    % sampling frequency
sens_in=1;      % input sensitivity [V/A]
sens_out=102*1e-3/9.81;     % output sensitivity [V/g]

% compute relevant parameters

n=max(size(y));    % number of samples
dt=1/fsamp;    % sampling time
t_end=(n-1)*dt;    % final time istant
t=0:dt:t_end;  % time vector

%% Data conversion

% convert from voltage to force and acceleration using sensitivities

x=x./sens_in;   % now x is in [A]
y=y./sens_out;  % now y is in [m/s^2]


%% Plot signals in time domain

if plot_time==1

figure

subplot(2,1,1)  % plot input force
plot(t,x)
grid on
title('Input force in time domain')
xlabel('t [s]')
ylabel('x [A]')
axis tight

subplot(2,1,2)  % plot output
plot(t,y)
grid on
title('Output in time domain')
xlabel('t [s]')
ylabel('y [m/s^2]')
axis tight

end

%% Signals in frequency domain

% compute the spectra of the signals

X=fft(x);
Y=fft(y);

% convert in only positive frequencies

[X_pos,f_pos]=positive_spectrum(X,fsamp);
Y_pos=positive_spectrum(Y,fsamp);

% plot input force in frequency domain

if plot_freq==1

figure

subplot(2,1,1)  % amplitude
semilogy(f_pos,abs(X_pos))
grid on
title('Input force in frequency domain')
xlabel('f [Hz]')
ylabel('|X|')
axis tight

subplot(2,1,2)  % phase
plot(f_pos,rad2deg(angle(X_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠X  [deg]')
axis tight

% plot output

figure

subplot(2,1,1)  % amplitude
semilogy(f_pos,abs(Y_pos))
grid on
title('Output in frequency domain')
xlabel('f [Hz]')
ylabel('|Y|')
axis tight

subplot(2,1,2)  % phase
plot(f_pos,rad2deg(angle(Y_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠Y  [deg]')
axis tight

end

%% T=10, 66% OVERLAP

% overlap and time duration of the subrecords used to average the power and
% cross spectra


T=10;       % assign time width of the subrecords
overlap=0.66;  % assign overlap of the subrecords

%% Power spectra and cross spectra (actually PSD)

% compute power and cross spectra

[Sxx,f_T]=crossSpectrum(x,x,T,fsamp,overlap);
Syy=crossSpectrum(y,y,T,fsamp,overlap); 
Sxy=crossSpectrum(x,y,T,fsamp,overlap);

% convert for only positive frequencies

[Gxx,f_T_pos]=positive_spectrum(Sxx,fsamp);
Gyy=positive_spectrum(Syy);
Gxy=positive_spectrum(Sxy);

% plot them

if plot_powcross==1

figure   % power spectrum of input
semilogy(f_T_pos,Gxx) 
xlabel('f [Hz]')
ylabel('G_x_x [A^2]')
grid on
title('Power spectrum of the input (T=10, ovrl=0.66)')
axis tight

figure  % power spectrum of output
semilogy(f_T_pos,Gyy)
xlabel('f [Hz]')
ylabel('G_y_y [(m/s^2)^2]')
grid on
title('Power spectrum of the output (T=10, ovrl=0.66)')
axis tight

figure  % cross spectrum

subplot(2,1,1)  % magnitude
semilogy(f_T_pos,abs(Gxy))
xlabel('f [Hz]')
ylabel('|G_x_y|')
grid on
title('Cross spectrum (T=10, ovrl=0.66)')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(Gxy)))
xlabel('f [Hz]')
ylabel('∠G_x_y  [deg]')
grid on
axis tight

end

% compute PSD: normalise for delta_f, since signals are random, so we need to
% take in account for leakage

delta_f=1/T;
Sxx=Sxx./delta_f;
Syy=Syy./delta_f;
Sxy=Sxy./delta_f;

% convert for only positive frequencies

[Gxx,f_T_pos]=positive_spectrum(Sxx,fsamp);
Gyy=positive_spectrum(Syy);
Gxy=positive_spectrum(Sxy);

% plot them

if plot_powcross==1

figure   % power spectrum of input
semilogy(f_T_pos,Gxx) 
xlabel('f [Hz]')
ylabel('G_x_x [A^2]')
grid on
title('PSD of the input (T=10, ovrl=0.66)')
axis tight

figure  % power spectrum of output
semilogy(f_T_pos,Gyy)
xlabel('f [Hz]')
ylabel('G_y_y [(m/s^2)^2]')
grid on
title('PSD of the output (T=10, ovrl=0.66)')
axis tight

figure  % cross spectrum

subplot(2,1,1)  % magnitude
semilogy(f_T_pos,abs(Gxy))
xlabel('f [Hz]')
ylabel('|G_x_y|')
grid on
title('Cross spectrum density (T=10, ovrl=0.66)')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(Gxy)))
xlabel('f [Hz]')
ylabel('∠G_x_y  [deg]')
grid on
axis tight

end

%% Coherence function

% compute coherence function

coh=coherence(Sxx,Syy,Sxy);

% plot it

if plot_coh==1

figure
plot(f_T_pos,coh);
grid on
xlabel('f [Hz]')
ylabel('\gamma^2_x_y')
title('Coherence function (T=10, ovrl=0.66)')
axis tight

end

%% Transfer function estimation

Gyx=conj(Gxy);
H1_pos=Gxy./Gxx;
H2_pos=Gyy./Gyx;

% plot H1 and H2

if plot_tf==1

figure  % H1

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H1_pos))
grid on
title('H_1 estimator (T=10, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('|H_1|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_1  [deg]')
axis tight

figure  % H2

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H2_pos))
grid on
title('H_2 estimator (T=10, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('|H_2|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_2  [deg]')
axis tight

end

%% Comparisons - H1, H2, and coherence

% plot H1, H2 and the coherence together

if plot_tf==1

figure

subplot(3,1,1)  % coherence
plot(f_T_pos,coh,'k')
grid on
title('Coherence (T=10, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('\gamma_x_y^2')
axis tight

subplot(3,1,2)  % H1 and H2 amplitude
semilogy(f_T_pos,abs(H1_pos))
hold on
semilogy(f_T_pos,abs(H2_pos))
grid on
title('Comparison between H_1 and H_2 (T=10, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('amplitude')
legend('H_1','H_2')
axis tight
hold off
    
subplot(3,1,3)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
hold on
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('phase [deg]')
legend('H_1','H_2')
axis tight
hold off

end

%% T=50, 66% OVERLAP
% overlap and time duration of the subrecords used to average the power and
% cross spectra


T=50;       % assign time width of the subrecords
overlap=0.66;  % assign overlap of the subrecords

%% Power spectra and cross spectra (actually PSD)

% compute power and cross spectra

[Sxx,f_T]=crossSpectrum(x,x,T,fsamp,overlap);
Syy=crossSpectrum(y,y,T,fsamp,overlap); 
Sxy=crossSpectrum(x,y,T,fsamp,overlap);

% convert for only positive frequencies

[Gxx,f_T_pos]=positive_spectrum(Sxx,fsamp);
Gyy=positive_spectrum(Syy);
Gxy=positive_spectrum(Sxy);

% plot them

if plot_powcross==1

figure   % power spectrum of input
semilogy(f_T_pos,Gxx) 
xlabel('f [Hz]')
ylabel('G_x_x [A^2]')
grid on
title('Power spectrum of the input (T=50, ovrl=0.66)')
axis tight

figure  % power spectrum of output
semilogy(f_T_pos,Gyy)
xlabel('f [Hz]')
ylabel('G_y_y [(m/s^2)^2]')
grid on
title('Power spectrum of the output (T=50, ovrl=0.66)')
axis tight

figure  % cross spectrum

subplot(2,1,1)  % magnitude
semilogy(f_T_pos,abs(Gxy))
xlabel('f [Hz]')
ylabel('|G_x_y|')
grid on
title('Cross spectrum (T=50, ovrl=0.66)')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(Gxy)))
xlabel('f [Hz]')
ylabel('∠G_x_y  [deg]')
grid on
axis tight

end

% compute PSD: normalise for delta_f, since signals are random, so we need to
% take in account for leakage

delta_f=1/T;
Sxx=Sxx./delta_f;
Syy=Syy./delta_f;
Sxy=Sxy./delta_f;

% convert for only positive frequencies

[Gxx,f_T_pos]=positive_spectrum(Sxx,fsamp);
Gyy=positive_spectrum(Syy);
Gxy=positive_spectrum(Sxy);

% plot them

if plot_powcross==1

figure   % power spectrum of input
semilogy(f_T_pos,Gxx) 
xlabel('f [Hz]')
ylabel('G_x_x [A^2]')
grid on
title('PSD of the input (T=50, ovrl=0.66)')
axis tight

figure  % power spectrum of output
semilogy(f_T_pos,Gyy)
xlabel('f [Hz]')
ylabel('G_y_y [(m/s^2)^2]')
grid on
title('PSD of the output (T=50, ovrl=0.66)')
axis tight

figure  % cross spectrum

subplot(2,1,1)  % magnitude
semilogy(f_T_pos,abs(Gxy))
xlabel('f [Hz]')
ylabel('|G_x_y|')
grid on
title('Cross spectrum density (T=50, ovrl=0.66)')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(Gxy)))
xlabel('f [Hz]')
ylabel('∠G_x_y  [deg]')
grid on
axis tight

end

%% Coherence function

% compute coherence function

coh=coherence(Sxx,Syy,Sxy);

% plot it

if plot_coh==1

figure
plot(f_T_pos,coh);
grid on
xlabel('f [Hz]')
ylabel('\gamma^2_x_y')
title('Coherence function (T=50, ovrl=0.66)')
axis tight

end

%% Transfer function estimation

Gyx=conj(Gxy);
H1_pos=Gxy./Gxx;
H2_pos=Gyy./Gyx;

% plot H1 and H2

if plot_tf==1

figure  % H1

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H1_pos))
grid on
title('H_1 estimator (T=50, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('|H_1|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_1  [deg]')
axis tight

figure  % H2

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H2_pos))
grid on
title('H_2 estimator (T=50, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('|H_2|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_2  [deg]')
axis tight

end

%% Comparisons - H1, H2, and coherence

% plot H1, H2 and the coherence together

if plot_tf==1

figure

subplot(3,1,1)  % coherence
plot(f_T_pos,coh,'k')
grid on
title('Coherence (T=50, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('\gamma_x_y^2')
axis tight

subplot(3,1,2)  % H1 and H2 amplitude
semilogy(f_T_pos,abs(H1_pos))
hold on
semilogy(f_T_pos,abs(H2_pos))
grid on
title('Comparison between H_1 and H_2 (T=50, ovrl=0.66)')
xlabel('f [Hz]')
ylabel('amplitude')
legend('H_1','H_2')
axis tight
hold off
    
subplot(3,1,3)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
hold on
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('phase [deg]')
legend('H_1','H_2')
axis tight
hold off

end