% ESE 3510 — Case Study 2: Fluorescence Lifetime Analysis
% Authors: Samir Afsary, Alex Ge, Ahmad Hamzeh

% CONTENTS:
%   Task 1a — Time-domain signals (excitation, fluorescence, detection)
%   Task 1b — Frequency response H_f(jw): numerical vs. theoretical
%   Task 1c — Effect of detector IRF on sampled signal y[n]
%   Task 2a — Why 100 kHz cannot recover 80 MHz fluorescence (aliasing)
%   Task 2b — Heterodyne detection: mixing to a recoverable beat frequency
%   Task 2c — Phase/modulation lifetime estimation from heterodyne signal
%   Task 3a — Phasor equations verified on single-exponential IRF
%   Task 3b — Per-pixel phasor analysis of real FLIM data
%

close all; clear; clc;

%%

% Parameters 
tau     = 1e-9;         % Fluorescence lifetime
T_L     = 12.5e-9;      % Laser period 
f_L     = 1 / T_L;      % Laser frequency 
omega_L = 2*pi * f_L;   % Laser angular frequency 
T_S     = 0.2e-9;       % ADC sampling period [ns]
dt      = 0.001e-9;      % small DT to simulate CT

%% Task 1a 

t = 0 : dt : 10*T_L - dt;   

% Create excitation impulse train
x = zeros(size(t));
for k = 0 : 9
    idx = round(k * T_L / dt) + 1;
    if idx <= length(t)
        x(idx) = 1 / dt;
    end
end

% Fluorescence IRF 
h_f = exp(-t / tau);     

% Calculate emitted fluorescence
f_full = conv(x, h_f) * dt;
f = f_full(1 : length(t));

% Detector IRF
sigma_t = 0.1e-9;
t0      = 0.175e-9;
h_d = exp(-(t - t0).^2 / sigma_t^2);

% Detector output
d_full = conv(f, h_d) * dt;
d = d_full(1 : length(t));

% Analog-to-Digital Conversion
step_size  = round(T_S / dt);
y_n        = d(1 : step_size : end);
t_discrete = t(1 : step_size : end);

% Plot functions
fig1 = figure('Name','Task 1a: Time-Domain Signals','NumberTitle','off', ...
              'Position',[50 50 900 750]);
% Excitation
subplot(5,1,1);
stem(t * 1e9, x * dt, 'filled', 'MarkerSize', 4, 'Color', [0.2 0.4 0.8]);
xlabel('Time (ns)'); ylabel('x(t)');
title('x(t): Excitation');
grid on; xlim([0,120]);

% Emission
subplot(5,1,2);
plot(t * 1e9, f, 'b', 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('f(t)');
title('f(t): Emitted Fluorescence');
grid on; xlim([0,120]);

% Continuous Detector Output
subplot(5,1,3);
plot(t * 1e9, d, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('d(t)');
title('d(t): Detector Output');
grid on; xlim([0,120]);

% Continuous Output
subplot(5,1,4);
plot(t * 1e9, d, 'Color', [0.6 0.1 0.6], 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('y(t)');
title('y(t): Continuous System Output');
grid on; xlim([0,120]);

% Digitized Output
subplot(5,1,5);
stem(t_discrete * 1e9, y_n, 'filled', 'MarkerSize', 4, 'Color', [0.6 0.1 0.6]);
xlabel('Time (ns)'); ylabel('y[n]');
title('y[n]: Digitized System Output');
grid on; xlim([0,120]);

sgtitle('Task 1a: Time-Domain Signals', 'FontSize', 13, 'FontWeight', 'bold');


%% Task 1b
dt_fine = dt / 100; 

% Recreate the time vector and the signal using the finer resolution
max_time = (length(h_f) - 1) * dt; 
t_fine = 0 : dt_fine : max_time; 
h_f_fine = exp(-t_fine/tau); 

% Numerical Frequency Response
N = length(h_f_fine);
H_f_num = fftshift(fft(h_f_fine)) * dt_fine; 
f_axis = (-N/2 : N/2-1) / (N * dt_fine);   

% Theoretical Frequency Response
omega_axis = 2*pi * f_axis;  
H_f_theory = tau ./ (1 + 1j * omega_axis * tau);

fig2 = figure('Name','Task 1b: Frequency Response H_f(jw)', ...
    'NumberTitle','off', 'Position',[50 50 900 550]);

% Magnitude
subplot(2,1,1);
plot(f_axis, abs(H_f_num),    'b-',  'LineWidth', 2.0); hold on;
plot(f_axis, abs(H_f_theory), 'r--', 'LineWidth', 1.5);
% Zoom in to ignore the edges of the fine simulation the limit is where
% the graph cuts of if you use dt instead of dt_fine
xlim([-5e11, 5e11]); 
xlabel('Frequency (Hz)'); 
ylabel('|H_f(j\omega)|');
title('Magnitude');
legend('Numerical FFT', 'Theoretical');
grid on;

% Phase
subplot(2,1,2);
plot(f_axis, angle(H_f_num)*180/pi,    'b-',  'LineWidth', 2.0); hold on;
plot(f_axis, angle(H_f_theory)*180/pi, 'r--', 'LineWidth', 1.5);
xlim([-5e11, 5e11]); % Zoom in to ignore the edges of the fine simulation
xlabel('Frequency (Hz)'); 
ylabel('Phase (degrees)');
title('Phase');
legend('Numerical FFT', 'Theoretical');
grid on;

sgtitle('Task 1b: Fluorescence Frequency Response H_f(j\omega)', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% Task 1c

%             fast,     medium,  slow   detectors
sigma_vals = [0.1e-9,   2.5e-9,  10.0e-9];
t0_vals =    [0.175e-9, 4e-9,    17e-9];
colors = {[0 0.4 1], [0.1 0.7 0.1], [0.8 0.1 0.1]};
labels = {'\sigma_t=0.1 ns, t_0=0.175 ns (fast)',  ...
          '\sigma_t=2.5 ns, t_0=4 ns (medium)', ...
          '\sigma_t=10.0 ns, t_0=17 ns (slow)'};

fig3 = figure('Name','Task 1c: Detector IRF Impact','NumberTitle','off', ...
              'Position',[50 50 1100 600]);

for k = 1 : 3
    % Build detector IRF
    h_dk = exp(-(t - t0_vals(k)).^2 / sigma_vals(k)^2);

    % Convolve fluorescence with this detector
    dk_full = conv(f, h_dk) * dt;
    dk = dk_full(1 : length(t));
    
    % Digitalize
    step_size  = round(T_S / dt);
    y_n        = dk(1 : step_size : end);
    t_discrete = t(1 : step_size : end);

    % Detector IRF
    subplot(2, 3, k);
    plot(t, h_dk, 'Color', colors{k}, 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('h_d(t)');
    title(labels{k}); grid on;
    if k == 1
        xlim([0,4e-10]);
    elseif k == 2
        xlim([0,1e-8]);
    else
        xlim([0,4e-8]);
    end

    % y[n]
    subplot(2, 3, k+3);
    stem(t_discrete, y_n, 'Color', colors{k});
    xlabel('Time (s)'); ylabel('y[n]');
    title('y[n]');
    grid on;
    xlim([0,125e-9]);
end

sgtitle({'Detector IRF Impact on y[n]'}, 'FontSize', 12, ...
    'FontWeight', 'bold');
