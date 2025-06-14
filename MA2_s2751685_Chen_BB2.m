% Phasing effect with a time-varying 1st-order All-Pass Filter chain.
% This implementation uses four All-Pass Filters (APFs), each associated 
% with a break frequency (notch frequency). 
% The notch frequencies vary over time, modulated by four independent 
% Low-Frequency Oscillators (LFOs).
% Each LFO is characterized by its own depth and rate, controlling the 
% extent and speed of frequency modulation.
% 
% Author: Zuyu Chen
% Date: 04/12/2024
clc
clear 
close all

Fs = 44100;                % Sample rate
D = [0.9 0.8 0.7 0.6]';    % LFO modulation depth
f = [0.2 0.35 0.25 0.2]';  % LFO rate

% Break frequencies 
fb_0 = 400;               % the first break/notch frequency 
fb = fb_0*[1 2 2^2 2^3]'; % exponentially spaced notch frequency series 

% fixed depth control
g = 1;

% Read in an input WAV file and store it in the vector 

% fn = 'Godin4_44.wav';
fn = 'Guitar_dry.wav';
% fn = 'KS_example.wav';
[x, ~] = audioread(['audio_samples/' fn]);

% Time varying notch frequencies
w_b = 2*pi*fb/Fs.*(1 + D.*sin(2*pi*f.*(0:length(x) - 1)/Fs));

% corresponding pole locations (approximated and normalized)
p_d = 1 - w_b; 

% Combine stereo to mono chanel
x = sum(x,2)/2;

% Preallocate the output vector y[n]
y = zeros(1, length(x));

% Preallocate the delay line buffer for 4 stages
z = zeros(4+1, length(x));  % Each row corresponds to a stage in the all-pass filter chain

for n = 2:length(x)
    % Assign the current input sample to the input of stage I
    z(1, n) = x(n);
    
    for i = 1:4
        % Compute the output of the current stage (i)
        % y(n-1): The output of the previous stage at the previous sample
        % x(n): The input of the current stage at the current sample
        % z(i, n): The input to the current stage at the current sample
        % z(i+1, n-1): The output of the previous stage at the last sample
        yStage = p_d(i, n) * z(i+1, n-1) + p_d(i, n) * z(i, n) - z(i, n-1);
        
        % Store the output of the current stage for use as input to the next stage
        z(i+1, n) = yStage;
    end
    
    % Combine the final stage's output with the input signal (scaled by gain g)
    y(n) = yStage + g * x(n);
end

soundsc(y,Fs) % After the effect
soundsc(x,Fs) % Before the effect

% Audio Full-Scale Normalization 
y = y - mean(y);
y = y/max(abs(y));
% Save the audio file
audiowrite('Phaser.wav', y, Fs);
