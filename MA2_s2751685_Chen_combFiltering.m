% This script will perform feedforward and feedback comb filtering 
% on a mono input WAV file x. 
% 
% Author: Zuyu Chen
% Date 29/11/2024

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
y_fb = zeros(length(x), 1); % feedback

% Preallocate the delay line buffer with M zeros
ff_dlinebuff = zeros(M, 1);     % feedforward delay line buffer
fb_dlinebuff = zeros(M, 1);     % feedback delay line buffer
% Comb filtering for both the feedforward and the feedback
for n = 1: length(x)
    
    % feedforward comb filtering
    y_ff(n) = x(n) + g*ff_dlinebuff(M);
    % update the feedforward delay line buffer
    ff_dlinebuff = circshift([ff_dlinebuff(1:end-1); x(n)], 1);
    % feedback comb filtering
    y_fb(n) = x(n) - g*fb_dlinebuff(M);
    % update the feedback delay line buffer
    fb_dlinebuff = circshift([fb_dlinebuff(1:end-1); y_fb(n)], 1);
     
end

soundsc(y_ff, Fs)
soundsc(y_fb, Fs)
soundsc(x, Fs)