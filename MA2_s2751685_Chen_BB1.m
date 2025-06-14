% This script implements a pointer-like circular buffer for feedforward
% combfiltering and flanger
% Author: Zuyu Chen
% Date 29/11/2024
%% ============Feedforward Comb Filter with Circular Buffer================

clc
clear
close all

% Set up the parameters
Fs = 44100;  % Sample rate
T = 0.025;   % Delay time in sec
g = 0.78;     % Comb filter coefficient (the 'strength' of the effect)
             % 
M = round(Fs*T);   % Delay in samples

% Read in an input WAV file and store it in the vector 

% fn = 'A440Hz.wav';
fn = 'Cath_cut.wav';
% fn = 'birchcanoe.wav';
% fn = 'bodhran-cutM.wav';
% fn = 'Godin4_44.wav';
% fn = 'KS_example.wav';
% fn = 'lathe.wav';
% fn = 'myks.wav';
[x, Fs] = audioread(['audio_samples/' fn]);

% Combine stereo to mono chanel
x = sum(x,2)/2;

% Preallocate the output vectors with the same length as the input vector
y_ff = zeros(length(x), 1); % feedforward 

% Preallocate the delay line buffer with M zeros
ff_dlinebuff = zeros(M, 1);     % feedforward delay line buffer

% Comb filtering for both the feedforward and the feedback
writeIdx = 1; % Write pointer

for n = 1: length(x)
    
   % Read from buffer (delayed sample)
    readIdx = writeIdx;
    delayedSample = ff_dlinebuff(readIdx);
    
    % Apply feedforward comb filter equation
    y_ff(n) = x(n) + g * delayedSample;
    
    % Write current input into buffer
    ff_dlinebuff(writeIdx) = x(n);
    
    % Increment write pointer (wrap around)
    writeIdx = mod(writeIdx, M) + 1;
     
end

soundsc(y_ff, Fs)
soundsc(x, Fs)

%% ======================Flanger with Circular Buffer======================
clc
clear
close all

% Set up the parameters

Fs = 44100;  % Sample rate
T = 0.025;   % Average delay time in sec
M0 = Fs*T;   % Average delay line length (n.d.)
g = 1;     % Comb filter coefficient (the 'strength' of the effect)
f0 = 0.2;    % LFO rate (Hz)
a = 1;       % LFO depth


% Read in an input WAV file and store it in the vector 

% fn = 'Cath_cut.wav';
fn = 'Guitar_dry.wav';
[x, Fs] = audioread(['audio_samples/' fn]);

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
writeIdx = 1;
for n = 1:length(x)

    % Read from buffer (delayed sample) 
    readIndx = round(mod(writeIdx - M(n) - 1, M0)) + 1;
    delayedSample = dlinebuff(readIndx);
    % Apply flanger equation
    y(n) = x(n) + g * delayedSample;
    % Write current input into buffer
    dlinebuff(writeIdx) = x(n);
    % Increment wirte pointer (wrap around)
    writeIdx = round(mod(writeIdx, M0)) + 1;

end
soundsc(y, Fs)
soundsc(x, Fs)

