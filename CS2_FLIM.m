%Fluorescence lifetime data from 2-photon excitation of NAD(P)H in breast cancer cells
load("FLIMhistogram.mat"); 
% x-y-t where x and y are spatial coordinates with 512 pixels each
% t is averaged over a few hundred laser periods with 80 MHz excitation
% time resolution and time vector given below
timeRes = 0.1280; % Time bin width in ns
time = timeRes*((1:size(FLIMhistogram,3))-1);

% fun thing to try if time - spatial binning


% estimate g and s values


% estimate lifetime - one value per x-y pixel


% visualize lifetime image 512x512


% create phasor plot - one point per pixel


% fun thing to try - can you segment out areas on the phasor plot and see
% where in the cells they correspond to?