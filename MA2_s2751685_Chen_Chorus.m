%% 
% This script will perform a chorus effect, based on feedforward 
% time-varying comb filtering of a mono input WAV file x 
% 
% Author: Zuyu Chen
% Date: 29/11/2024

clc
clear
close all

% Params
% mode = 'round';
mode = 'interpolation';
N = 4;
Q = 100;

Fs0 = 44100; % Original sample rate 
% Fs1 = 44800; % Changed sample rate
Fs1 = 22050; % Changed sample rate
T1 = 0.009;  % Average delay time of buffer 1 in sec
PO1 = Fs1*T1;% Average delay line length of buffer 1 in samples
g1 = 1;      % Effect strength of buffer 1
f1 = 1;      % LFO1 rate (Hz)
t1 = 0.0023; % LFO1 depth in time
D1 = t1*Fs1; % LFO1 depth in samples

T2 = 0.001;  % Average delay time of buffer 2 in sec
PO2 = Fs1*T2; % Average delay line length of buffer 2 in samples
g2 = 1;    % Effect strength of buffer 2
f2 = 10;     % LFO2 rate (Hz)
t2 = 0.00025; % LFO2 depth in time
D2 = t2*Fs1;  % LFO2 depth in samples

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
M1 = PO1 + D1*sin(2*pi*f1*(0:length(x) - 1)/Fs1);
M2 = PO2 + D2*sin(2*pi*f2*(0:length(x) - 1)/Fs1);

% Preallocate the output vector y[n]
y = zeros(1, length(x));

% Preallocate the delay line buffer with M+D zeros
tolerance = 1e-10;  % apply a small number to prevent rounding to zeros
dlinebuff1 = zeros(ceil(PO1+D1+tolerance), 1);
dlinebuff2 = zeros(ceil(PO2+D2+tolerance), 1);
if isequal(mode,'interpolation')
    dlinebuff1 = [dlinebuff1; zeros(N,1)];
    dlinebuff2 = [dlinebuff2; zeros(N,1)];
end

% perform the sample-wise filtering operation 
for n = 1:length(x)
    
    switch mode
        case 'round'
            % just round the delay line index (no interpolation)
            y(n) = x(n) + g1 * dlinebuff1(ceil(M1(n)+tolerance))...
                        + g2 * dlinebuff2(ceil(M2(n)+tolerance));
        case 'interpolation'
            
            % Apply Lagrange Interpolation
            % =======================================
        
            T = MA2_s2751685_Chen_Linterp(N, Q, 1);
            a = (-Q/2 + (0:Q-1))/Q ;  % sample locations
            
            % round the delay line index to the nearest position in the Lagrange
            % interpolation samples
            a_M1 = M1(n) - floor(M1(n)) - 1/2;   % -1/2 <= a_M1 < 1/2
            [~, row_idx] = min(abs(a_M1 - a));

            % 
            val1 = T(row_idx, :) * dlinebuff1(floor(M1(n)) - (N/2-1): ...
                                              floor(M1(n)) - (N/2-1) + N -1);
            
            a_M2 = M2(n) - floor(M2(n)) - 1/2;   % -1/2 <= a_M2 < 1/2
            [~, row_idx] = min(abs(a_M2 - a));
            val2 = T(row_idx, :) * dlinebuff2(floor(M2(n)) - (N/2-1): ...
                                              floor(M2(n)) - (N/2-1) + N -1);
            
            % calculate the output
            y(n) = x(n) + g1*val1 + g2*val2;
    end
   
    
    % update the delay line buffers
    dlinebuff1 = [x(n); dlinebuff1(1:end-1)];
    dlinebuff2 = [x(n); dlinebuff2(1:end-1)];

end
% Resample to convert 
y = resample(y, Fs1, Fs0);
x = resample(x, Fs1, Fs0);
soundsc(y, Fs1)
soundsc(x, Fs1)
% Audio Full-Scale Normalization 
y = y - mean(y);
y = y/max(abs(y));

% audiowrite('Linterp_N=10,Q=100.wav', y, Fs1);
% audiowrite('Round.wav', y,  Fs1)

