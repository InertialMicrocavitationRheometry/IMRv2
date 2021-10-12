function [] = f_fourier(Tout,Rout,Ro_w,lm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Fs = length(Rout)/Tout(end);    % Sampling frequency
T = 1/Fs;                       % Sample time
L = length(Rout);               % Length of signal
t = (0:L-1)*T;                  % Time vector

Rc = (Rout-mean(Rout))/Ro_w-1;

NFFT = floor(2^(nextpow2(L)+5)); % Next power of 2 from length of y
Y = fft(Rc,NFFT)/L;
fvec = Fs/2*linspace(0,1,NFFT/2+1);
% Plot single-sided amplitude spectrum
plot(fvec,2*abs(Y(1:NFFT/2+1)),lm,'LineWidth',2) 
end