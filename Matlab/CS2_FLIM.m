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

%% Task 1
% *Part a* 

t = 0 : T_S : 5*T_L - T_S;   

% Create excitation impulse train
x = zeros(size(t));
for k = 0 : floor(5*T_L / T_L) - 1
    idx = round(k * T_L / T_S) + 1;
    if idx <= length(t)
        x(idx) = 1 / T_S;
    end
end

% Fluorescence IRF 
h_f = exp(-t / tau);     

% Calculate emitted fluorescence
f_full = conv(x, h_f) * T_S;
f = f_full(1 : length(t));

% Detector IRF
sigma_t = 0.3e-9;   % Detector temporal width [ns]  (fast PMT ~300 ps)
t0 = 1e-9;     % Detector time delay [ns]
h_d = zeros(size(t));
mask = (t >= t0);
h_d(mask) = exp(-(t(mask) - t0) / sigma_t);   % Zero before t0 (causal)

% Detector output
d_full = conv(f, h_d) * T_S;
d = d_full(1 : length(t));

% Sampled signal y[n]: t is already on a T_S grid, so y[n] = d(t) 
y = d;

% Plot functions
fig1 = figure('Name','Task 1a: Time-Domain Signals','NumberTitle','off', ...
              'Position',[50 50 900 750]);

subplot(5,1,1);
stem(t, x * T_S, 'filled', 'MarkerSize', 5, 'Color', [0.2 0.4 0.8]);
xlabel('Time (s)'); ylabel('x(t)');
title('Excitation x(t): Impulse Train at 80 MHz'); grid on;

subplot(5,1,2);
plot(t, f, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('f(t)');
title(sprintf('Emitted Fluorescence f(t): Single-Exp Decay, \\tau = %.1f ns', tau));
grid on;

subplot(5,1,3);
plot(t, h_d, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('h_d(t)');
title(sprintf('Detector IRF h_d(t): \\sigma_t = %.1f ns, t_0 = %.1f ns', sigma_t, t0));
grid on;

subplot(5,1,4);
plot(t, d, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('d(t)');
title('Detector Output d(t) = f(t) * h_d(t)  (broadened by detector)'); grid on;

subplot(5,1,5);
stem(t, y, 'filled', 'MarkerSize', 2, 'Color', [0.6 0.1 0.6]);
xlabel('Time (s)'); ylabel('y[n]');
title(sprintf('Sampled Signal y[n] at T_S = %.1f ns  (5 GS/s ADC)', T_S)); grid on;

sgtitle('Task 1a: Time-Domain Signal Chain', 'FontSize', 13, 'FontWeight', 'bold');

%% 
% *Part b*

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

%% 
% *Part C*

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

%% Task 2
% *Part a*
f_L_hz     = 80e6;                       % 80 MHz [Hz]
omega_L_hz = 2*pi * f_L_hz;             % [rad/s]
tau_s      = tau * 1e-9;                 % Convert ns -> s for Hz-domain math

% Modulation depth and phase of the fluorescence
M_theory   = 1 / sqrt(1 + (omega_L_hz * tau_s)^2);
phi_theory = -atan(omega_L_hz * tau_s);  % [rad]

fprintf('\n=== Task 2a ===\n');
fprintf('Laser frequency       f_L  = %g MHz\n', f_L_hz/1e6);
fprintf('Modulation depth      M    = %.4f\n', M_theory);
fprintf('Phase shift           phi  = %.4f rad = %.2f deg\n', ...
        phi_theory, phi_theory*180/pi);
fprintf('100 kHz Nyquist             = 50 kHz << 80 MHz  => ALIASED\n');
fprintf('Conclusion: slow sampling loses phase => cannot estimate tau\n\n');

% Show aliasing visually over a 0.5 µs window
T_S_slow  = 1 / 100e3;                  % 10 µs per sample
t_fast_s  = 0 : 1e-10 : 5e-7;          % 0.5 µs at 10 GHz (true signal)
t_slow_s  = 0 : T_S_slow : 5e-7;       % 100 kHz samples

x_fast    = 1 + cos(omega_L_hz * t_fast_s);
f_fast    = 1 + M_theory * cos(omega_L_hz * t_fast_s + phi_theory);
x_slow_s  = 1 + cos(omega_L_hz * t_slow_s);
f_slow_s  = 1 + M_theory * cos(omega_L_hz * t_slow_s + phi_theory);

fig4 = figure('Name','Task 2a: Aliasing at 100 kHz','NumberTitle','off', ...
              'Position',[50 50 900 500]);

subplot(2,1,1);
plot(t_fast_s*1e6, x_fast, 'b', 'LineWidth', 1); hold on;
stem(t_slow_s*1e6, x_slow_s, 'r', 'filled', 'MarkerSize', 7);
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Excitation x(t): True 80 MHz vs. 100 kHz Samples (aliased)');
legend('True x(t)', '100 kHz samples'); grid on; xlim([0, 0.5]);

subplot(2,1,2);
plot(t_fast_s*1e6, f_fast, 'b', 'LineWidth', 1); hold on;
stem(t_slow_s*1e6, f_slow_s, 'r', 'filled', 'MarkerSize', 7);
xlabel('Time (\mus)'); ylabel('Amplitude');
title('Fluorescence f(t): 80 MHz aliases at 100 kHz => phase info destroyed');
legend('True f(t)', '100 kHz samples'); grid on; xlim([0, 0.5]);

sgtitle('Task 2a: 100 kHz Sampling Cannot Recover 80 MHz Phase', ...
        'FontSize', 12, 'FontWeight', 'bold');

%  TASK 2a: Why 100 kHz Cannot Recover Fluorescence Lifetime
%
%  Excitation: x(t) = 1 + cos(omega_L * t),  f_L = 80 MHz
%  Fluorescence for single exponential (no IRFs):
%    f(t) = 1 + M*cos(omega_L*t + phi)
%  where:
%    M   = 1 / sqrt(1 + (omega_L * tau)^2)   [modulation depth]
%    phi = -atan(omega_L * tau)               [phase shift, radians]
%
%  Nyquist for 100 kHz ADC = 50 kHz << 80 MHz => severe aliasing.
%  The 80 MHz signal folds to an alias whose frequency depends on the
%  ratio 80 MHz / 100 kHz = 800, giving NO usable phase information.
% =========================================================================


%% =========================================================================
%  TASK 2b: Heterodyne Detection
%
%  Multiply f(t) by g(t) = cos(omega_g * t) where omega_g is slightly
%  offset from omega_L so a low-frequency "beat" falls within the
%  100 kHz ADC bandwidth.
%
%  Mixing product (using trig identity):
%    f(t)*g(t) ~ M/2 * cos((omega_L - omega_g)*t + phi) + high-freq terms
%  where the beat frequency  delta_f = f_L - f_g  is << 50 kHz (Nyquist).
%
%  After a low-pass filter (LPF), only the beat survives.
%  The beat PRESERVES the phase phi (and modulation M) of the fluorescence,
%  so lifetime can be recovered via heterodyne comparison with the
%  beat from a reference copy of the excitation.
% =========================================================================

delta_f   = 1e3;                         % Intermediate/beat frequency [Hz] = 1 kHz
f_g_hz    = f_L_hz - delta_f;           % Mixer (heterodyne) frequency [Hz]
omega_g   = 2*pi * f_g_hz;

% Simulate using 1 µs time steps over 5 ms (many beat cycles).
% 1 µs >> 1/80 MHz so the 80 MHz carrier is not resolved — that is fine;
% we only need to track the envelope (beat), which is at 1 kHz.
T_S_mid   = 1e-6;                        % 1 µs intermediate sampling [s]
t_mid     = 0 : T_S_mid : 5e-3;         % 5 ms window

f_mid     = 1 + M_theory * cos(omega_L_hz * t_mid + phi_theory);  % Fluorescence
x_mid     = 1 + cos(omega_L_hz * t_mid);                          % Excitation
g_mid     = cos(omega_g * t_mid);                                   % Mixer

det_f     = f_mid .* g_mid;     % Fluorescence after mixing
det_x     = x_mid .* g_mid;     % Excitation after mixing (reference channel)

% LPF: zero out FFT bins above 2 kHz to isolate the beat
N_mid     = length(t_mid);
f_mid_ax  = (-N_mid/2 : N_mid/2-1) / (N_mid * T_S_mid);   % Centered freq axis [Hz]
LPF_cut   = 2 * delta_f;   % 2 kHz cutoff

F_det_f   = fftshift(fft(det_f));
F_det_x   = fftshift(fft(det_x));
lpf_mask  = (abs(f_mid_ax) <= LPF_cut);

y_beat_f  = real(ifft(ifftshift(F_det_f .* lpf_mask)));
y_beat_x  = real(ifft(ifftshift(F_det_x .* lpf_mask)));

% ---- Figure 5: Task 2b — Time Domain ----
fig5 = figure('Name','Task 2b: Heterodyne - Time Domain','NumberTitle','off', ...
              'Position',[50 50 900 650]);

subplot(3,1,1);
plot(t_mid*1e3, f_mid, 'b', 'LineWidth', 0.6);
xlabel('Time (ms)'); ylabel('f(t)'); grid on;
title(sprintf('Fluorescence f(t) at %g MHz — oscillates too fast to sample', f_L_hz/1e6));

subplot(3,1,2);
plot(t_mid*1e3, det_f, 'Color',[0.1 0.6 0.1], 'LineWidth', 0.6);
xlabel('Time (ms)'); ylabel('f(t) \cdot g(t)'); grid on;
title(sprintf('After Mixing: f(t)\\cdotg(t), g at f_g = %g MHz', f_g_hz/1e6));

subplot(3,1,3);
plot(t_mid*1e3, y_beat_f, 'r',  'LineWidth', 1.8); hold on;
plot(t_mid*1e3, y_beat_x, 'b--','LineWidth', 1.8);
xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
title(sprintf('After LPF: Beat at \\Delta f = %g kHz — detectable at 100 kHz', delta_f/1e3));
legend('Fluorescence beat','Excitation beat (ref)');

sgtitle('Task 2b: Heterodyne Detection — Time Domain', ...
        'FontSize', 12, 'FontWeight', 'bold');

% ---- Figure 6: Task 2b — Frequency Domain ----
fig6 = figure('Name','Task 2b: Heterodyne - Freq Domain','NumberTitle','off', ...
              'Position',[50 50 900 650]);

subplot(3,1,1);
F_f_mid = fftshift(fft(f_mid));
plot(f_mid_ax/1e6, abs(F_f_mid)/N_mid, 'b', 'LineWidth', 1.2);
xlabel('Frequency (MHz)'); ylabel('|F(\omega)|');
title('Spectrum of f(t): peaks at \pm80 MHz');
xlim([-200, 200]); grid on;

subplot(3,1,2);
plot(f_mid_ax/1e6, abs(F_det_f)/N_mid, 'Color',[0.1 0.6 0.1], 'LineWidth', 1.2);
xlabel('Frequency (MHz)'); ylabel('|F(\omega)|');
title(sprintf('After Mixing: peaks at f_L \\pm f_g = \\pm%g MHz and 0 MHz', ...
              (f_L_hz + f_g_hz)/2/1e6));
xlim([-200, 200]); grid on;

subplot(3,1,3);
plot(f_mid_ax/1e3, abs(F_det_f .* lpf_mask)/N_mid, 'r',  'LineWidth', 1.8); hold on;
plot(f_mid_ax/1e3, abs(F_det_x .* lpf_mask)/N_mid, 'b--','LineWidth', 1.8);
xlabel('Frequency (kHz)'); ylabel('|Y(j\omega)|');
title(sprintf('After LPF: Beat at \\pm%g kHz — within 100 kHz bandwidth', delta_f/1e3));
legend('|Y_f(j\omega)|','|Y_x(j\omega)| (ref)');
xlim([-5, 5]); grid on;

sgtitle('Task 2b: Heterodyne Detection — Frequency Domain', ...
        'FontSize', 12, 'FontWeight', 'bold');


%% =========================================================================
%  TASK 2c: Phase Shift and Modulation — Heterodyne Lifetime Estimation
%
%  Phase relationship:
%    Phase of fluorescence beat = phi_f = phi_fluorescence + phi_g
%    Phase of excitation beat   = phi_x = 0              + phi_g
%    => Delta_phi = phi_f - phi_x = phi_fluorescence = -atan(omega_L * tau)
%    => tau = -tan(Delta_phi) / omega_L
%
%  Modulation relationship:
%    Delta_M = M_fluorescence / M_excitation = 1 / sqrt(1 + (omega_L*tau)^2)
%    => tau = sqrt(1/Delta_M^2 - 1) / omega_L
% =========================================================================

% Resample beats at 100 kHz (well above the 1 kHz beat frequency)
T_S_100k  = 1 / 100e3;
t_100k    = t_mid(1) : T_S_100k : t_mid(end);

y_f_100k  = interp1(t_mid, y_beat_f, t_100k, 'linear', 0);
y_x_100k  = interp1(t_mid, y_beat_x, t_100k, 'linear', 0);

% Extract phase and amplitude at the beat frequency using DFT
N_100k   = length(t_100k);
Y_f_dft  = fft(y_f_100k);
Y_x_dft  = fft(y_x_100k);

% DFT bin index for delta_f
k_beat   = round(delta_f * N_100k * T_S_100k) + 1;

phase_f    = angle(Y_f_dft(k_beat));
phase_x    = angle(Y_x_dft(k_beat));
delta_phi  = phase_f - phase_x;    % Phase difference [rad]

amp_f      = abs(Y_f_dft(k_beat));
amp_x      = abs(Y_x_dft(k_beat));
delta_M    = amp_f / amp_x;        % Modulation ratio

% Lifetime estimates
tau_from_phase = -tan(delta_phi) / omega_L_hz * 1e9;   % [ns]
tau_from_mod   = sqrt(max(0, 1/delta_M^2 - 1)) / omega_L_hz * 1e9;  % [ns]

fprintf('=== Task 2c ===\n');
fprintf('Phase difference   Delta_phi = %.4f rad = %.2f deg\n', ...
        delta_phi, delta_phi*180/pi);
fprintf('Modulation ratio   Delta_M   = %.4f\n', delta_M);
fprintf('Lifetime (phase):  tau = %.4f ns  (true: %.1f ns)\n', tau_from_phase, tau);
fprintf('Lifetime (modulation): tau = %.4f ns  (true: %.1f ns)\n\n', tau_from_mod, tau);

% ---- Figure 7: Task 2c ----
fig7 = figure('Name','Task 2c: Phase Shift and Lifetime Estimation', ...
              'NumberTitle','off','Position',[50 50 900 500]);

subplot(2,1,1);
plot(t_100k*1e3, y_x_100k, 'b', 'LineWidth', 1.8); hold on;
plot(t_100k*1e3, y_f_100k, 'r', 'LineWidth', 1.8);
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Beat Signals at 100 kHz — Phase Shift \\Delta\\phi = %.2f°', ...
              delta_phi*180/pi));
legend('Excitation beat (reference)','Fluorescence beat'); grid on;
xlim([0, 5]);

subplot(2,1,2);
f_ax_100k = (0 : N_100k-1) / (N_100k * T_S_100k);
plot(f_ax_100k/1e3, abs(Y_x_dft)/N_100k, 'b', 'LineWidth', 1.8); hold on;
plot(f_ax_100k/1e3, abs(Y_f_dft)/N_100k, 'r', 'LineWidth', 1.8);
xlabel('Frequency (kHz)'); ylabel('Magnitude');
title(sprintf('Spectra — \\tau from phase: %.3f ns, from modulation: %.3f ns', ...
              tau_from_phase, tau_from_mod));
legend('Excitation','Fluorescence'); grid on;
xlim([0, 5]);

sgtitle('Task 2c: Heterodyne Lifetime Estimation', ...
        'FontSize', 12, 'FontWeight', 'bold');


%% =========================================================================
%  TASK 3a: Phasor Equations — Verification on Single Exponential
%
%  Phasor components (integrated over one laser period T_L):
%    g = integral( h_f(t)*cos(omega_L*t) dt ) / integral( h_f(t) dt )
%    s = integral( h_f(t)*sin(omega_L*t) dt ) / integral( h_f(t) dt )
%
%  Analytical solution for h_f(t) = exp(-t/tau) integrated over [0, inf]:
%    g = 1 / ( 1 + (omega_L*tau)^2 )
%    s = (omega_L*tau) / ( 1 + (omega_L*tau)^2 )
%
%  Lifetime recovery:
%    tau = s / (omega_L * g) = (omega_L*tau) / (omega_L * 1) = tau  [checks out]
%
%  Interpretation:
%    g = "in-phase" component — higher g means shorter lifetime (less phase lag)
%    s = "quadrature" component — peaks at omega_L*tau = 1  (tau = 1/(2*pi*f_L))
%    Together (g,s) lie on the universal semicircle of radius 0.5 centered at (0.5,0)
% =========================================================================

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


%% =========================================================================
%  TASK 3b: Phasor Analysis of Real FLIM Data
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