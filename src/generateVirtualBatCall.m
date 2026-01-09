function call = generateVirtualBatCall(f0, f1, d, fs, tail)
% Generates a virtual bat call, simulating a frequency-modulated (FM) sweep
% commonly used by bats for echolocation. The call is designed to start at
% a frequency of f1, linearly sweep down to a frequency of fo over the duration
% of d, and is sampled at a rate of fs. The call is windowed with a
% Hanning window to smooth start and end transients, and an optional
% spectral emphasis can be applied at the midpoint frequency to simulate
% bat call characteristics.
% 
% Inputs:
%   f0      - Start frequency of the call in Hz
%   f1      - End frequency of the call in Hz
%   d       - Duration of the call in milliseconds
%   fs      - Sampling rate in Hz
%   tail    - Zero padding as % of d, 0 for no padding. Adds before and
%             after the signal
%
% Output:
%   call    - The generated virtual bat call signal, including pre- and
%             post-silence padding.
%
% Example:
%   call = generateVirtualBatCall(f0, f1, d, fs);
%   call = generateVirtualBatCall(25000, 30000, 0.002, 192000);
%   This generates a 2 ms bat call sweeping from 30 kHz to 25 kHz,
%   sampled at 192 kHz.

% Calculate the frequency for spectral emphasis (fmax)
% The fmax calculation aims to emphasize a frequency in the call's
% spectrum. The formula provided is an example; modify as needed
% to match specific bat call characteristics.
fmax = mean([f0,f1])-f0/3;

% Create the linear chirp signal
% t defines the time vector for the duration of the call
t = 0:1/fs:d-1/fs;
% Generate a linear chirp with frequencies sweeping from f0 to f1
vel = chirp(t, f0, d, f1, 'quadratic');
% Reverse the chirp to create a descending frequency sweep
vel = fliplr(vel);

% Generate the filter for spectral emphasis at fmax
% Define the frequency bands and magnitudes for the yulewalk filter
fb = [0 2*[f0 fmax f1] fs]./fs; % Normalized frequency bands
m = [0 0 1 0 0];                % Magnitude response at specified bands
% Generate the yulewalk filter coefficients
[yb, ya] = yulewalk(4, fb, m);
% Create an impulse response of the filter
[h, ~] = impz(yb, ya, fs/1000);
fmaxir = h;

% Apply a Hanning window to the sweep signal to smooth transitions
window = hanning(length(vel));
vel = vel.*window';

% Filter the signal with the designed filter to emphasize fmax
% It is important to filter after windowing to avoid spectral leakage
% vel = filtfilt(fmaxir, 1, vel);
% vel = vel./max(abs(vel));

% Add silence before and after the call to simulate natural pauses
call = [zeros(round((d*fs)*tail/100), 1); vel'; zeros(round((d*fs)*tail/100), 1)];
end
