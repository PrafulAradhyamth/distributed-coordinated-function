% Combined Wi-Fi Saturation Throughput Model (Simulation & Analytical)
%
% This script calculates and plots the saturation throughput of an IEEE 802.11
% Distributed Coordination Function (DCF) network using two methods:
% 1. A simplified discrete-event simulation.
% 2. Bianchi's analytical model.
%
% This allows for direct comparison of the simulated performance against
% the theoretical upper bound.

clear; close all; clc; % Clear workspace, close figures, clear command window

%% 1. Define Common System Parameters (from Bianchi's Parameter Table)
% These parameters are common to both the simulation and the analytical model.

payload = 8184;      % Payload size in bits (e.g., 1024 bytes * 8 bits/byte)
MACheader = 272;     % MAC header size in bits (e.g., 34 bytes * 8 bits/byte)
PHYheader = 128;     % PHY header size in bits (e.g., 16 bytes * 8 bits/byte)
ack = 112;           % ACK frame size in bits (excluding PHY header)
ProDelay = 1;        % Propagation Delay in microseconds (µs) - common in analytical models

% Derived packet sizes
packet_data = payload + MACheader + PHYheader; % Total data packet size (bits)
packet_ack = ack + PHYheader;                  % Total ACK packet size (bits)

% Interframe Spaces and Slot Time (µs)
SIFS = 28;           % Short Interframe Space (µs)
SlotTime = 50;       % Slot Time (µs)
DIFS = 128;          % DCF Interframe Space (µs)

% Derived time durations in Slot Times for analytical model
% Ts: Time for a successful transmission
% This includes: Data Packet Tx + SIFS + ACK Tx + DIFS
Ts = (packet_data + SIFS + ProDelay + packet_ack + DIFS + ProDelay) / SlotTime;

% Tc: Time for a collision
% This includes: Data Packet Tx (collided) + DIFS
Tc = (packet_data + DIFS + ProDelay) / SlotTime;

%% 2. Define Contention Window (CW) and Maximum Backoff Stage Scenarios
% W: Initial Contention Window (CWmin + 1, as Bianchi's W is 0-indexed)
% m: Maximum backoff stage (CW doubles 'm' times)

W_scenarios = [32, 32, 128]; % Initial Contention Window values (W_0)
m_scenarios = [3, 5, 3];     % Max backoff stages

% Number of stations to simulate/analyze
num_stations_range = 3:1:50;
num_data_points = length(num_stations_range);

%% 3. Simulation and Analytical Model Loop
% Iterate through each scenario (different W and m values)

figure; % Create a new figure for the plot
hold on; % Keep the plot open for multiple lines
grid on;
xlabel('Number of Stations');
ylabel('Saturation Throughput (Mbps)');
title('Saturation Throughput: Simulation vs. Analytical Model');

legend_entries = {}; % To store legend strings

for j = 1 : length(W_scenarios)
    current_W = W_scenarios(j);
    current_m = m_scenarios(j);

    % Arrays to store results for the current scenario
    throughput_sim = zeros(1, num_data_points);
    throughput_analytical = zeros(1, num_data_points);

    fprintf('--- Running Scenario: W=%d, m=%d ---\n', current_W, current_m);

    for n_idx = 1:num_data_points
        n = num_stations_range(n_idx); % Current number of stations

        %% Discrete-Event Simulation (first MATLAB code)
        simutime_sim = 0; % Total simulation time for this run
        Spkt_sim = 0;     % Number of successfully transmitted packets

        % Initialize station states and waiting times
        StationStage = zeros(1, n); % Current backoff stage for each station
        waitingtime = zeros(1, n);  % Current backoff timer for each station (in µs)

        for i = 1 : n
            StationStage(i) = 0; % Start at backoff stage 0
            % Initial backoff: DIFS + random value within CW_0
            waitingtime(i) = DIFS + floor(current_W * rand) * SlotTime;
        end

        % Main simulation loop: run until 10,000 successful packets
        while (Spkt_sim < 10000)
            % Find the minimum waiting time (when the next event occurs)
            idle_end = min(waitingtime);

            % Count transmissions occurring at this time
            num_tran = sum(waitingtime == idle_end);

            % Update simulation time by the idle period
            simutime_sim = simutime_sim + idle_end;

            if num_tran == 1 % Only one transmission: SUCCESSFUL
                Spkt_sim = Spkt_sim + 1;

                % Add time for successful packet exchange
                simutime_sim = simutime_sim + packet_data + SIFS + packet_ack + DIFS;

                for i = 1 : n
                    if waitingtime(i) == idle_end
                        % Successful station: reset to stage 0, new random CW
                        StationStage(i) = 0;
                        waitingtime(i) = floor(current_W * rand) * SlotTime;
                    else
                        % Other stations: decrement their timers by idle_end
                        waitingtime(i) = waitingtime(i) - idle_end;
                    end
                end
            else % num_tran > 1: COLLISION
                % Add time for collision (full packet transmission) and DIFS
                simutime_sim = simutime_sim + packet_data + DIFS;

                for i = 1 : n
                    if waitingtime(i) == idle_end % Colliding stations
                        % Increase backoff stage if not at max
                        if StationStage(i) < current_m
                            StationStage(i) = StationStage(i) + 1;
                        end
                        % New random backoff from doubled CW
                        CW_current = (2^StationStage(i)) * current_W;
                        waitingtime(i) = floor(CW_current * rand) * SlotTime;
                    else
                        % Other stations: decrement their timers by idle_end
                        waitingtime(i) = waitingtime(i) - idle_end;
                    end
                end
            end
        end
        % Calculate throughput for simulation (Mbps)
        throughput_sim(n_idx) = Spkt_sim * payload / simutime_sim / 1e6; % bits/µs -> Mbps


        %% Bianchi's Analytical Model (second MATLAB code)

        % Define the function to find 'p' (collision probability)
        % Note: Bianchi's W is CW_min+1. Here, current_W is already CW_min+1.
        fun_p = @(p) (p - (1 - ( (2 * (1 - 2*p)) / ((1 - 2*p) * (current_W + 1) + p * current_W * (1 - (2*p)^current_m)) ) )^(n - 1) );

        % Find the root of the equation for p (collision probability)
        % Use a numerical solver (fzero) within the interval [0, 1)
        % A robust way to find the root, ensuring 0 <= p < 1
        % Check if a root exists in a small range to avoid errors when fzero fails
        options = optimset('Display','off'); % Suppress fzero output
        try
            pj = fzero(fun_p, [0, 0.99999999], options); % Search in [0, 1)
        catch
            % If fzero fails to find a root (e.g., for very low n),
            % it implies very low collision probability.
            % A common workaround is to assume p is very small or zero.
            % For this specific Bianchi function, a root often exists but might
            % be outside the default search range or require a different initial guess.
            % Let's try a different approach if fzero fails on the initial range.
            % For robustness, we might need a more sophisticated root-finding strategy.
            % For now, set to a small value if it fails to avoid breaking the script.
            warning('fzero could not find a root for p. Setting p to 0.001 for n=%d', n);
            pj = 0.001; % Fallback
        end


        % Calculate tae (probability a station transmits in a given slot)
        tae = (2 * (1 - 2*pj)) / ((1 - 2*pj) * (current_W + 1) + pj * current_W * (1 - (2*pj)^current_m));
        
        % In some edge cases for small N or large W, tae can become complex or
        % negative if p is not precisely found or the formula's domain limits are hit.
        % Ensure tae is a real and non-negative probability.
        if ~isreal(tae) || tae < 0
             tae = 0; % Or handle as an error/special case
        end
        
        % Ptr: Probability of at least one transmission in the slot
        Ptr = 1 - (1 - tae)^n;

        % Ps: Probability of a successful transmission, given there's a transmission
        if Ptr == 0
            Ps = 0; % No transmissions, no success
        else
            Ps = (n * tae * (1 - tae)^(n - 1)) / Ptr;
        end

        % Efi: Expected number of idle slots before a transmission
        if Ptr == 0
            Efi = 0; % No transmissions, no idle before transmission
        else
            Efi = (1 / Ptr) - 1;
        end

        % Calculate throughput for analytical model (Mbps)
        % Ensure Efi is non-negative and the denominator is not zero
        denominator = Efi + Ps * Ts + (1 - Ps) * Tc;
        if denominator <= 0
            throughput_analytical(n_idx) = 0; % Avoid division by zero or negative throughput
        else
            throughput_analytical(n_idx) = (Ps * payload / SlotTime) / denominator / 1e6; % bits/µs -> Mbps
        end
        
        fprintf('  n=%d: Sim Thr=%.2f Mbps, Analytical Thr=%.2f Mbps\n', n, throughput_sim(n_idx), throughput_analytical(n_idx));

    end

    % Plot results for the current scenario
    plot(num_stations_range, throughput_sim, '-', 'LineWidth', 2);
    legend_entries{end+1} = sprintf('Sim: W=%d, m=%d', current_W, current_m);

    plot(num_stations_range, throughput_analytical, '--', 'LineWidth', 1.5);
    legend_entries{end+1} = sprintf('Analytical: W=%d, m=%d', current_W, current_m);

end

% Finalize plot
legend(legend_entries, 'Location', 'best');
hold off;

fprintf('\nSimulation and Analytical Model complete. Check the plot.\n');