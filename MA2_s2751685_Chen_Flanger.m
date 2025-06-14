%% 
% This script will perform a flanger effect, based on feedforward 
% time-varying comb filtering of a mono input WAV file x 
% 
% Author: Zuyu Chen
% Date: 29/11/2024

clc
clear
close all

% Set up the parameters

Fs = 44100;  % Sample rate
T = 0.025;   % Average delay time in sec
M0 = Fs*T;   % Average delay line length (n.d.)
g = 1;       % Comb filter coefficient (the 'strength' of the effect)
f0 = 0.2;    % LFO rate (Hz)
a = 1;       % LFO depth


% Read in an input WAV file and store it in the vector 

% fn = 'A440Hz.wav';
% fn = 'Cath_cut.wav';
% fn = 'birchcanoe.wav';
% fn = 'bodhran-cutM.wav';
% fn = 'Godin4_44.wav';
fn = 'Guitar_dry.wav';
% fn = 'KS_example.wav';
% fn = 'lathe.wav';
% fn = 'myks.wav';
[x, ~] = audioread(['audio_samples/' fn]);

% Combine stereo to mono chanel
x = sum(x,2)/2;

% Pre-compute the delay profile vector M[n]
M = M0*(1 + a*sin(2*pi*f0*(0:length(x) - 1)/Fs));

% Preallocate the output vector y[n]
y = zeros(1, length(x));

% Preallocate the delay line buffer with 2*M0 zeros
tolerance = 1e-10;  % apply a small number to prevent rounding to zeros
dlinebuff = zeros(ceil(2*M0+tolerance), 1);

% perform the filtering operation sample-wise
for n = 1:length(x)
    
    y(n) = x(n) + g * dlinebuff(ceil(M(n)+tolerance)); 
    dlinebuff = [x(n); dlinebuff(1:end-1)];

end
soundsc(y, Fs)
soundsc(x, Fs)

% Audio Full-Scale Normalization 
% y = y - mean(y);
% y = y/max(abs(y));
% audiowrite('flanger_round.wav', y, Fs);
