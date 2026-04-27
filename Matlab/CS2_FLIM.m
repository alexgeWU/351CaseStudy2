% ESE 3510 — Case Study 2: Fluorescence Lifetime Analysis
% Authors: Samir Afsary, Alex Ge, Ahmad Hamzeh

% CONTENTS:
%   Task 1a — Time-domain signals (excitation, fluorescence, detection)
%   Task 1b — Frequency response H_f(jw): numerical vs. theoretical
%   Task 1c — Effect of detector IRF on sampled signal y[n]
%   Task 2a — Why 100 kHz cannot recover 100 MHz fluorescence (aliasing)
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

%% Task 1a - Time-domain signals (excitation, fluorescence, detection)

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
sigma_t = 0.3e-9;   % Detector temporal width [ns]  (fast PMT ~300 ps)
t0      = 1e-9;     % Detector time delay [ns]
h_d     = zeros(size(t));
mask    = (t >= t0);
h_d(mask) = exp(-(t(mask) - t0).^2 / sigma_t^2);   % Zero before t0 (causal)

% Detector output
d_full = conv(f, h_d) * dt;
d = d_full(1 : length(t));

% Analog-to-Digital Conversion
step_size  = round(T_S / dt);
y_n        = d(1 : step_size : end);
t_discrete = t(1 : step_size : end); % Time vector for plotting the discrete points

% Plot functions
fig1 = figure('Name','Task 1a: Time-Domain Signals','NumberTitle','off', ...
              'Position',[50 50 900 750]);
% Excitation
subplot(5,1,1);
stem(t * 1e9, x * dt, 'filled', 'MarkerSize', 4, 'Color', [0.2 0.4 0.8]);
xlabel('Time (ns)'); ylabel('x(t)');
title('Excitation x(t): Impulse Train at 80 MHz');
grid on; xlim([0,120]);

% Emission
subplot(5,1,2);
plot(t * 1e9, f, 'b', 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('f(t)');
title(sprintf('Emitted Fluorescence f(t): Single-Exp Decay, \\tau = %.1f ns', tau * 1e9));
grid on; xlim([0,120]);

% Continuous Detector Output
subplot(5,1,3);
plot(t * 1e9, d, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('d(t)');
title('Continuous Detector Output d(t) = f(t) * h_d(t)');
grid on; xlim([0,120]);

% Continuous Output
subplot(5,1,4);
plot(t * 1e9, d, 'Color', [0.6 0.1 0.6], 'LineWidth', 1.5);
xlabel('Time (ns)'); ylabel('y[n]');
title('Continuous Output y(t) = d(t)');
grid on; xlim([0,120]);

% Digitized Output
subplot(5,1,5);
% Overlay the continuous signal lightly behind the discrete samples
stem(t_discrete * 1e9, y_n, 'filled', 'MarkerSize', 4, 'Color', [0.6 0.1 0.6]);
xlabel('Time (ns)'); ylabel('y[n]');
title(sprintf('Sampled Signal y[n] at T_S = %.1f ns  (5 GS/s ADC)', T_S * 1e9));
grid on; xlim([0,120]);

sgtitle('Task 1a: Time-Domain Signal Chain', 'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(gcf, '../docs/figs/1atimedomain.jpg');

%% Task 1b - Frequency response H_f(jw): numerical vs. theoretical

% Numerical frequency response via DFT
N = length(h_f);
H_f_num = fft(h_f) * T_S; % Scale by T_S for continuous-FT approximation
f_axis = (0 : N-1) / (N * T_S);   

% Theoretical frequency response at the same frequencies
omega_axis = 2*pi * f_axis;  
H_f_theory = tau ./ (1 + 1j * omega_axis * tau);

% Plot frequency responses
fig2 = figure('Name','Task 1b: Frequency Response H_f(jw)', ...
    'NumberTitle','off', 'Position',[50 50 900 550]);

subplot(2,1,1);
plot(f_axis, abs(H_f_num),    'b-',  'LineWidth', 2.0); hold on;
plot(f_axis, abs(H_f_theory), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (GHz)'); ylabel('|H_f(j\omega)|');
title('Magnitude — Numerical FFT vs. Theoretical \tau/(1+j\omega\tau)');
legend('Numerical FFT', 'Theoretical');
grid on;

subplot(2,1,2);
plot(f_axis, angle(H_f_num)*180/pi,    'b-',  'LineWidth', 2.0); hold on;
plot(f_axis, angle(H_f_theory)*180/pi, 'r--', 'LineWidth', 1.5);
xlabel('Frequency (GHz)'); ylabel('Phase (degrees)');
title('Phase — Numerical FFT vs. Theoretical  \angle H_f = -atan(\omega\tau)');
legend('Numerical FFT', 'Theoretical');
grid on;

sgtitle('Task 1b: Fluorescence Frequency Response H_f(j\omega)', ...
        'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(gcf, '../docs/figs/1bfreqresponse.jpg');

%% Task 1c - Effect of detector IRF on sampled signal y[n]

% Three detector scenarios: fast, medium, slow
sigma_vals = [0.1e-9,  2.5e-9,  10.0e-9]; % Detector widths
t0_vals = [0.3e-9,  1.0e-9,  2.5e-9];     % Time delays 
colors = {[0 0.4 1], [0.1 0.7 0.1], [0.8 0.1 0.1]};
labels = {'\sigma_t=0.1 ns (fast)',  ...
          '\sigma_t=2.5 ns (medium)', ...
          '\sigma_t=10.0 ns (slow)'};

fig3 = figure('Name','Task 1c: Detector IRF Impact','NumberTitle','off', ...
              'Position',[50 50 1100 600]);

for k = 1 : 3
    % Build detector IRF for this scenario
    h_dk = zeros(size(t));
    mask_k = (t >= t0_vals(k));
    h_dk(mask_k) = exp(-(t(mask_k) - t0_vals(k)) / sigma_vals(k));

    % Convolve fluorescence with this detector
    dk_full = conv(f, h_dk) * T_S;
    dk = dk_full(1 : length(t));

    % Top row: detector IRFs
    subplot(2, 3, k);
    plot(t, h_dk, 'Color', colors{k}, 'LineWidth', 1.8);
    xlabel('Time (s)'); ylabel('h_d(t)');
    title(labels{k}); grid on;

    % Bottom row: resulting y[n]
    subplot(2, 3, k+3);
    stem(t, dk, 'filled', 'MarkerSize', 2, 'Color', colors{k});
    xlabel('Time (s)'); ylabel('y[n]');
    title(sprintf('y[n]: \\sigma_t=%.1f ns, t_0=%.1f ns', ...
        sigma_vals(k) * 10^9, t0_vals(k) * 10^9));
    grid on;
end

sgtitle({'Detector IRF Impact on y[n]'}, 'FontSize', 12, ...
    'FontWeight', 'bold');
exportgraphics(gcf, '../docs/figs/1cdetectorchanges.jpg');

%% Task 2a - Why 100 kHz cannot recover 100 MHz fluorescence (aliasing)
f_L_hz = 100e6;               % Laser modulation frequency
omega_L_hz = 2*pi * f_L_hz;      

% Show aliasing visually over a 0.5 µs window
T_S_slow = 1e-5;                      
T_S_fast = 1e-11;
t_fast_s = 0 : T_S_fast : 2e-5;         
t_slow_s = 0 : T_S_slow : 2e-5;       

% Compute the output of the ground truth using convolution
x_fast = 1 + cos(omega_L_hz * t_fast_s);
x_slow = 1 + cos(omega_L_hz * t_slow_s);
h_f_fast = 1 / tau * exp(-t_fast_s / tau);
f_fast = conv(x_fast, h_f_fast) * T_S_fast;
f_fast = f_fast(1 : length(t_fast_s));

% Sample slower frequency emission from ground truth
step = round(T_S_slow / T_S_fast);
f_slow = f_fast(1:step:end);

% Calculate theoretical modulation depth and phase of the fluorescence
M = 1 / sqrt(1 + (omega_L_hz * tau)^2);
phi = -atan(omega_L_hz * tau); 

% Create upper/lower bounds for visualization
upper_bound = (1 + M) * ones(1, numel(t_fast_s));
lower_bound = (1 - M) * ones(1, numel(t_fast_s));

fig4 = figure('Name','Aliasing at 100 kHz','NumberTitle','off', ...
              'Position',[50 50 900 500]);

subplot(2,1,1);
plot(t_fast_s*1e6, x_fast, 'b', 'LineWidth', 1); hold on;
stem(t_slow_s*1e6, x_slow, 'r', 'filled', 'MarkerSize', 7);
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Excitation x(t): Ground Truth 100 MHz vs. Sampled 100 kHz');
legend('True x(t)', '100 kHz samples'); grid on; xlim([0, 0.5]);

subplot(2,1,2);
plot(t_fast_s*1e6, f_fast, 'b', 'LineWidth', 1); hold on;
stem(t_slow_s*1e6, f_slow, 'r', 'filled', 'MarkerSize', 7);
plot(t_fast_s*1e6, upper_bound, 'LineWidth', 2, 'Color', 'green');
plot(t_fast_s*1e6, lower_bound, 'LineWidth', 2, 'Color', 'green');
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Emission f(t)');
legend('True f(t)', '100 kHz samples', 'Amplitude Bounds'); grid on; 
xlim([0, 20]);

sgtitle(['100 kHz Sampling vs. Ground Truth of 100 MHz Excitation ' ...
    'and Emission Signals'], 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(gcf, '../docs/figs/2abound.jpg')

%% Task 2b - Heterodyne detection: mixing to a recoverable beat frequency

delta_f = 1e3;                         
f_n_hz = f_L_hz - delta_f;           
omega_n = 2*pi * f_n_hz;

T_S_mid = 1e-6;                      
t_mid = 0 : T_S_mid : 5e-3 - T_S_mid;         

f_mid = 1 + M * cos(omega_L_hz * t_mid + phi);  % Fluorescence
x_mid = 1 + cos(omega_L_hz * t_mid);            % Excitation
g = cos(omega_n * t_mid);                       % Mixer

mixed_f = (M/2)*cos((omega_L_hz - omega_n)*t_mid + phi);
mixed_x = (1/2)*cos((omega_L_hz - omega_n)*t_mid);

% FFT setup
N  = length(t_mid);
fs = 1/T_S_mid;
freq = (-N/2:N/2-1)*(fs/N);

Y_f = fftshift(fft(mixed_f));
Y_x = fftshift(fft(mixed_x));

% Time-domain plots
figure;
subplot(3,1,1);
plot(t_mid*1e3, g, 'LineWidth', 1.8);
xlabel('Time (ms)');
ylabel('g(t)');
title('Mixer Signal');
grid on;

subplot(3,1,2);
plot(t_mid*1e3, mixed_f, 'b', 'LineWidth', 1.8);
xlabel('Time (ms)');
ylabel('f_{det}(t)g(t)');
title('Low-Frequency Term After Mixing + LPF');
grid on;

subplot(3,1,3);
plot(t_mid*1e3, mixed_f, 'b', 'LineWidth', 1.8); hold on;
plot(t_mid*1e3, mixed_x, 'r', 'LineWidth', 1.8);
xlabel('Time (ms)');
ylabel('Amplitude');
title('After LPF: 1 kHz Beat Signal');
legend('Fluorescence beat', 'Reference beat');
grid on;

exportgraphics(gcf, '../docs/figs/2btimedomain.jpg')

% Frequency-domain plots
figure;
subplot(2,1,1);
plot(freq/1e3, abs(Y_f)/N, 'b', 'LineWidth', 1.8); hold on;
plot(freq/1e3, abs(Y_x)/N, 'r', 'LineWidth', 1.8);
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('Magnitude of Beat Signals');
legend('Fluorescence beat', 'Reference beat');
xlim([-5 5]);
grid on;

subplot(2,1,2);
plot(freq/1e3, angle(Y_f), 'b', 'LineWidth', 1.8); hold on;
plot(freq/1e3, angle(Y_x), 'r', 'LineWidth', 1.8);
xlabel('Frequency (kHz)');
ylabel('Phase [rad]');
title('Phase of Beat Signals');
legend('Fluorescence beat', 'Reference beat');
xlim([-5 5]);
grid on;

exportgraphics(gcf, '../docs/figs/2bfreqdomain.jpg')

%% Task 2c - Phase/modulation lifetime estimation from heterodyne signal

[~, k_beat] = min(abs(freq - delta_f));

phase_f = angle(Y_f(k_beat));
phase_x = angle(Y_x(k_beat));

delta_phi = angle(exp(1j*(phase_f - phase_x)));   % wrapped phase difference
delta_M   = abs(Y_f(k_beat)) / abs(Y_x(k_beat));  % modulation ratio

% Lifetime estimates
tau_from_phase = -tan(delta_phi) / omega_L_hz * 1e9;   
tau_from_mod   = sqrt(max(0, 1/delta_M^2 - 1)) / omega_L_hz * 1e9;

figure;
plot(t_mid*1e3, mixed_x, 'r', 'LineWidth', 1.8); hold on;
plot(t_mid*1e3, mixed_f, 'b', 'LineWidth', 1.8);
xlabel('Time (ms)');
ylabel('Amplitude');
title(sprintf('Beat Signals: \\Delta\\phi = %.2f^\\circ, \\Delta M = %.3f', ...
      delta_phi*180/pi, delta_M));
legend('Excitation beat (reference)', 'Fluorescence beat');
grid on;
xlim([0 5]);

exportgraphics(gcf, '../docs/figs/2clifetime.jpg')
%% TASK 3a: Phasor Equations — Verification on Single Exponential

% Analytical phasor for tau = 1 ns, omega_L in rad/ns
g_theory_3a = 1 / (1 + (omega_L * tau)^2);
s_theory_3a = (omega_L * tau) / (1 + (omega_L * tau)^2);
tau_from_phasor_theory = s_theory_3a / (omega_L * g_theory_3a);

% Numerical phasor: integrate h_f over one period using trapezoidal rule
t_phasor  = 0 : 0.001 : T_L;             % Fine grid [ns]
h_f_p     = exp(-t_phasor / tau);         % IRF samples

norm_hf   = trapz(t_phasor, h_f_p);
g_num_3a  = trapz(t_phasor, h_f_p .* cos(omega_L * t_phasor)) / norm_hf;
s_num_3a  = trapz(t_phasor, h_f_p .* sin(omega_L * t_phasor)) / norm_hf;
tau_from_phasor_num = s_num_3a / (omega_L * g_num_3a);

fprintf('=== Task 3a: Phasor Verification (tau = %.1f ns) ===\n', tau);
fprintf('Theoretical: g=%.4f  s=%.4f  tau=%.4f ns\n', ...
        g_theory_3a, s_theory_3a, tau_from_phasor_theory);
fprintf('Numerical:   g=%.4f  s=%.4f  tau=%.4f ns\n\n', ...
        g_num_3a, s_num_3a, tau_from_phasor_num);

% Universal semicircle: all single-exponential lifetimes lie on this curve
theta_sc = linspace(0, pi, 500);
g_sc = 0.5 + 0.5*cos(theta_sc);
s_sc = 0.5*sin(theta_sc);

% Mark several tau values on the semicircle for reference
tau_marks = [0.2, 0.4, 0.8, 1.0, 2.0, 4.0, 8.0];

fig8 = figure('Name','Task 3a: Phasor Diagram','NumberTitle','off', ...
              'Position',[50 50 700 500]);
plot(g_sc, s_sc, 'k-', 'LineWidth', 2); hold on;

for tm = tau_marks
    g_tm = 1 / (1 + (omega_L*tm)^2);
    s_tm = (omega_L*tm) / (1 + (omega_L*tm)^2);
    plot(g_tm, s_tm, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [0.7 0.7 0.7]);
    text(g_tm + 0.02, s_tm + 0.01, sprintf('%.1f ns', tm), 'FontSize', 8);
end

plot(g_theory_3a, s_theory_3a, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot(g_num_3a,    s_num_3a,    'b+', 'MarkerSize', 12, 'LineWidth', 2.5);

xlabel('g  (in-phase component)', 'FontSize', 11);
ylabel('s  (quadrature component)', 'FontSize', 11);
title(sprintf('Phasor Diagram — \\tau = %.1f ns lies on the universal semicircle', tau), ...
      'FontSize', 11);
legend('Universal semicircle','Lifetime labels','Theoretical','Numerical', ...
       'Location','northwest');
grid on; axis equal; xlim([0, 1]); ylim([0, 0.6]);

sgtitle('Task 3a: Phasor Components g and s', ...
        'FontSize', 13, 'FontWeight', 'bold');


%% TASK 3b: Phasor Analysis of Real FLIM Data

fprintf('=== Task 3b: Loading FLIM data ===\n');

% Try current directory first, then ../input/ for submission folder structure
try
    load('data\FLIMhistogram.mat');
    fprintf('Loaded FLIMhistogram from current directory.\n');
catch
    try
        load('../input/FLIMhistogram.mat');
        fprintf('Loaded FLIMhistogram from ../input/\n');
    catch
        error(['FLIMhistogram.mat not found. ', ...
               'Place it in the same folder as main.m and re-run.']);
    end
end

% Time axis (given in starter code)
timeRes_FLIM = 0.1280;    % Time bin width [ns]
time_FLIM    = timeRes_FLIM * (0 : size(FLIMhistogram,3)-1);  % [ns]

% Laser period and angular frequency for phasor computation
T_L_FLIM     = 12.5;                    % [ns]
omega_L_FLIM = 2*pi / T_L_FLIM;        % [rad/ns]

[Nx, Ny, Nt] = size(FLIMhistogram);
fprintf('Image size: %d x %d pixels, %d time bins (%.3f ns/bin)\n', ...
        Nx, Ny, Nt, timeRes_FLIM);

% Reshape to (Npix × Nt) matrix for vectorized phasor computation
H = reshape(double(FLIMhistogram), Nx*Ny, Nt);

% Precompute cosine and sine vectors at the fundamental laser frequency
cos_vec = cos(omega_L_FLIM * time_FLIM);   % 1×Nt
sin_vec = sin(omega_L_FLIM * time_FLIM);   % 1×Nt

% Total photon count per pixel (the normalization denominator)
norm_H = sum(H, 2);    % Npix×1

% Only compute phasor for pixels with enough photons to be meaningful
valid = (norm_H > 10);

g_vec = zeros(Nx*Ny, 1);
s_vec = zeros(Nx*Ny, 1);

% Vectorized dot products: each row of H times cos_vec/sin_vec
g_vec(valid) = (H(valid, :) * cos_vec') ./ norm_H(valid);
s_vec(valid) = (H(valid, :) * sin_vec') ./ norm_H(valid);

% Lifetime from phasor: tau = s / (omega_L * g)
tau_vec = s_vec ./ (omega_L_FLIM * g_vec);

% Reshape phasor components and lifetime back to image grids
g_map   = reshape(g_vec, Nx, Ny);
s_map   = reshape(s_vec, Nx, Ny);
tau_map = reshape(tau_vec, Nx, Ny);

% Mask pixels with unphysical lifetimes (empty, negative, or > 6 ns)
tau_map(~reshape(valid, Nx, Ny)) = NaN;
tau_map(tau_map < 0)  = NaN;
tau_map(tau_map > 6)  = NaN;

fprintf('Mean   lifetime (valid pixels): %.3f ns\n', mean(tau_map(:), 'omitnan'));
fprintf('Median lifetime (valid pixels): %.3f ns\n', median(tau_map(:), 'omitnan'));

% ---- Figure 9: Lifetime Image ----
fig9 = figure('Name','Task 3b: Fluorescence Lifetime Image','NumberTitle','off', ...
              'Position',[50 50 700 650]);

imagesc(tau_map, [0, 3]);   % Color range 0–3 ns (typical for NAD(P)H)
colormap(jet);
cb = colorbar;
cb.Label.String = 'Lifetime \tau (ns)';
cb.FontSize = 11;
axis image; axis off;
title({'Fluorescence Lifetime Image', ...
       'NAD(P)H in Breast Cancer Cells (2-photon FLIM)'}, 'FontSize', 12);

% ---- Figure 10: Sample Pixel Decay Curves ----
% Find the pixel with the highest photon count as a representative example
[~, idx_best] = max(norm_H);
[px_b, py_b]  = ind2sub([Nx, Ny], idx_best);

% Find a second pixel: nearest to the median lifetime among bright pixels
tau_valid_vec = tau_vec;
tau_valid_vec(~valid | tau_vec < 0 | tau_vec > 6) = NaN;
tau_median_val = median(tau_valid_vec(:), 'omitnan');
bright_idx     = find(valid & norm_H > 50);
[~, offset]    = min(abs(tau_valid_vec(bright_idx) - tau_median_val));
idx_med        = bright_idx(offset);
[px_m, py_m]   = ind2sub([Nx, Ny], idx_med);

fig10 = figure('Name','Task 3b: Sample Decay Curves','NumberTitle','off', ...
               'Position',[50 50 900 420]);

for kk = 1:2
    if kk == 1
        px = px_b; py = py_b; lbl = 'Brightest pixel';
    else
        px = px_m; py = py_m; lbl = 'Median-lifetime pixel';
    end
    subplot(1,2,kk);
    decay = squeeze(FLIMhistogram(px, py, :));
    plot(time_FLIM, decay, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
    tau_px = tau_map(px, py);
    if ~isnan(tau_px) && tau_px > 0
        A_fit = max(decay);
        plot(time_FLIM, A_fit * exp(-time_FLIM / tau_px), 'r--', 'LineWidth', 1.8);
        legend('Data', sprintf('Phasor \\tau = %.2f ns', tau_px), ...
               'Location','northeast');
    end
    xlabel('Time (ns)'); ylabel('Photon Count');
    title(sprintf('%s (%d,%d)', lbl, px, py)); grid on;
end

sgtitle('Task 3b: Sample Pixel Fluorescence Decay Curves', ...
        'FontSize', 12, 'FontWeight', 'bold');

% ---- Figure 11: Phasor Plot (density image) ----
% Extract valid phasor points in physical range
valid_pts = valid & g_vec > 0 & g_vec < 1 & s_vec >= 0 & s_vec < 0.6;
g_pts     = g_vec(valid_pts);
s_pts     = s_vec(valid_pts);

% Build 2D density histogram (log scale avoids bright-pixel dominance)
g_edges  = linspace(0, 1,   201);
s_edges  = linspace(0, 0.55, 151);
[counts, ~, ~] = histcounts2(g_pts, s_pts, g_edges, s_edges);
g_ctrs   = (g_edges(1:end-1) + g_edges(2:end)) / 2;
s_ctrs   = (s_edges(1:end-1) + s_edges(2:end)) / 2;

fig11 = figure('Name','Task 3b: Phasor Plot','NumberTitle','off', ...
               'Position',[50 50 800 580]);

imagesc(g_ctrs, s_ctrs, log1p(counts'));   % Transpose: x=g, y=s
set(gca, 'YDir', 'normal');
colormap(hot);
cb2 = colorbar;
cb2.Label.String = 'log(1 + pixel count)';
hold on;

% Overlay universal semicircle
plot(g_sc, s_sc, 'w-', 'LineWidth', 2);

% Mark lifetime values on the semicircle
tau_labels_flim = [0.2, 0.4, 0.8, 1.0, 1.6, 2.0, 3.2, 4.0];
for tl = tau_labels_flim
    g_tl = 1 / (1 + (omega_L_FLIM * tl)^2);
    s_tl = (omega_L_FLIM * tl) / (1 + (omega_L_FLIM * tl)^2);
    plot(g_tl, s_tl, 'wo', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
    text(g_tl + 0.01, s_tl + 0.01, sprintf('%.1f ns', tl), ...
         'Color', 'w', 'FontSize', 8);
end

xlabel('g  (in-phase)', 'FontSize', 11);
ylabel('s  (quadrature)', 'FontSize', 11);
title({'Task 3b: Phasor Plot — All Pixels', ...
       'Pixel density (hot colormap, log scale); white curve = universal semicircle'}, ...
      'FontSize', 11);
xlim([0, 1]); ylim([0, 0.55]); grid on;

sgtitle('Task 3b: FLIM Phasor Plot', 'FontSize', 13, 'FontWeight', 'bold');

fprintf('\nAll tasks complete — %d figures generated.\n', 11);
%
%  Data: FLIMhistogram.mat  —  512×512 spatial pixels, ~98 time bins
%  Each pixel contains a time-resolved photon count histogram h[n].
%  Laser: 80 MHz two-photon excitation => T_L = 12.5 ns.
%
%  Per-pixel phasor:
%    g_px = sum( h[n]*cos(omega_L*t[n]) ) / sum( h[n] )
%    s_px = sum( h[n]*sin(omega_L*t[n]) ) / sum( h[n] )
%    tau_px = s_px / (omega_L * g_px)
%
%  Expected results for NAD(P)H in cancer cells:
%    Free  NAD(P)H: tau ~ 0.3–0.5 ns
%    Bound NAD(P)H: tau ~ 1–4 ns
%    (Points cluster inside the semicircle for mixed/multi-exponential pixels)
% =========================================================================