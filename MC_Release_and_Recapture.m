%% Simulating Release and Recapture in Optical Tweezers using Monte Carlo
clear
clc
% Generating random seed numbers
rng('shuffle')
% Defining constants, all the units are in SI
m_Rb = 87e-3/6.02e23;
k_B = 1.38e-23;
g = 9.8;
lambda = 852e-9; % lambda of the trapping laser.

% Defining trapping parameters
trap_freq_trans = 2 * pi * 390e3;
trap_freq_longi = 2 * pi * 68e3;
trap_depth = 600e-6 * k_B;
omega = [trap_freq_longi, trap_freq_trans, trap_freq_trans];
waist = 1.1e-6;

% Defining scanning parameters
T_range = linspace(10e-6, 100e-6, 10); % Temperature range
% T = 30e-6;
release_times = linspace(0.1e-6, 100e-6, 1000); % Release time range
% release_time = 30e-6;

N_trials = 1000;
results = zeros(length(T_range), length(release_times), N_trials);
% Loop over Temperature
for i_T = 1:length(T_range)
    T = T_range(i_T);
    % Loop over release_time
    for i_release_time = 1:length(release_times)
        release_time = release_times(i_release_time);
        % Loop over exp_time
        for i_trial = 1:N_trials
            % Initialize position and velocity randomly
            % For the 3 directions, the first one is longitudinal (z), the other two are transverse (x, y)
            r0 = randn(1, 3) .* sqrt(k_B * T ./ (m_Rb * omega.^2));
            v0 = randn(1, 3) .* sqrt(k_B * T ./ m_Rb);
        
            % Calculating final position and kinetic energy
            r = r0 + v0 * release_time - [1/2, 0, 0] * g * release_time^2;
            v = v0 - [1, 0, 0] * g * release_time;
            kinetic_energy = 1/2 * m_Rb * sum(v.^2);
            U = calc_potential(trap_depth, waist, r, lambda);
            if kinetic_energy > U
                results(i_T, i_release_time, i_trial) = 0;
%                 disp("Lossed")
            else
                results(i_T, i_release_time, i_trial) = 1;
%                 disp("Recaptured")
            end
        end
    end
end
%% Result processing
recapture_probability = squeeze(sum(results, 3))/N_trials;
%% Visualize
figure
hold on
for i_T = 1:length(T_range)
    plot(release_times * 1e6, recapture_probability(i_T,:), 'LineWidth', 2, 'DisplayName', ['T = ', num2str(T_range(i_T) * 1e6), ' \muK'])
end
legend('show')
xlabel('Release time (\mus)')
ylabel('Recapture probability')
title('Recapture Probability Versus Release time for Different Temperatures')
hold off
set(gca, "FontSize", 18)

%% Funs
function U = calc_potential(trap_depth, waist, r, lambda)
    Z_R = pi * waist^2 / lambda;
    w = waist * sqrt(1 + (r(1)/Z_R)^2 ); % waist radius at z = r(1)
    % Calculating local light intensity ratio with the bottom of the trap
    ratio =  exp(-2*(r(2)^2 + r(3)^2)/w^2) / (pi * w^2)  * (pi * waist^2);
%     disp("ratio="+num2str(ratio))
%     disp([Z_R, w, waist, r])
    U = trap_depth * ratio;
end