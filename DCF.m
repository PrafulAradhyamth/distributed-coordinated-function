% Bianchi's Model: Simulation and Analytical Comparison in One Plot

% Define parameters
payload = 8184; 
MACheader = 272; 
PHYheader = 128; 
ack = 112;

pkt = payload + MACheader + PHYheader;
ACK = 112 + PHYheader;

SIFS = 28;
SlotTime = 50;
DIFS = 128;
ProDelay = 1;

packet = payload + MACheader + PHYheader;
Ts = (PHYheader + MACheader + payload + SIFS + ProDelay + ACK + DIFS + ProDelay) / SlotTime;
Tc = (PHYheader + MACheader + payload + DIFS + ProDelay) / SlotTime;

W = [32, 32, 128];
m = [3, 5, 3];

figure;
hold on;

colors = ['r', 'g', 'b']; % Different colors for each (W, m)

for j = 1:3
    sim_throughput = zeros(1, 48);
    ana_throughput = zeros(1, 48);

    for n = 3:50
        % ---- Simulation ----
        Spkt = 0;
        simutime = 0;
        StationStage = zeros(1, n);
        waitingtime = DIFS + floor(W(j) * rand(1, n)) * SlotTime;

        while Spkt < 10000
            idle_end = min(waitingtime);
            num_tran = sum(waitingtime == idle_end);

            if num_tran == 1
                Spkt = Spkt + 1;
                simutime = simutime + pkt + idle_end + SIFS + ACK + DIFS;

                for i = 1:n
                    if waitingtime(i) == idle_end
                        StationStage(i) = 0;
                        waitingtime(i) = floor((W(j) - 1) * rand) * SlotTime;
                    else
                        waitingtime(i) = waitingtime(i) - idle_end;
                    end
                end
            else
                simutime = simutime + idle_end + pkt + DIFS;

                for i = 1:n
                    if waitingtime(i) == idle_end
                        StationStage(i) = min(StationStage(i) + 1, m(j));
                        CW = 2^StationStage(i) * W(j);
                        waitingtime(i) = floor((CW - 1) * rand) * SlotTime;
                    else
                        waitingtime(i) = waitingtime(i) - idle_end;
                    end
                end
            end
        end

        sim_throughput(n - 2) = Spkt * payload / simutime;

        % ---- Analytical ----
        fun = @(p) p - 1 + (1 - (2*(1 - 2*p)) / ((1 - 2*p)*(W(j) + 1) + p * W(j) * (1 - (2*p)^m(j))))^(n - 1);
        pj = fzero(fun, [0, 1]);
        tau = (2 * (1 - 2 * pj)) / ((1 - 2 * pj)*(W(j) + 1) + pj * W(j) * (1 - (2 * pj)^m(j)));
        Ptr = 1 - (1 - tau)^n;
        Ps = (n * tau * (1 - tau)^(n - 1)) / Ptr;
        Efi = 1 / Ptr - 1;

        ana_throughput(n - 2) = (Ps * payload / SlotTime) / (Efi + Ps * Ts + (1 - Ps) * Tc);
    end

    % Plot both simulation and analytical results
    n_range = 3:50;
    plot(n_range, sim_throughput, ['-' colors(j)], 'LineWidth', 2);         % solid line for simulation
    plot(n_range, ana_throughput, ['--' colors(j)], 'LineWidth', 1.5);      % dashed line for analytical
end

% Plot formatting
title('Saturation Throughput: Simulation vs Analytical');
xlabel('Number of Stations');
ylabel('Throughput (bits per time unit)');
legend({ ...
    'Sim: W=32, m=3', 'Analytical: W=32, m=3', ...
    'Sim: W=32, m=5', 'Analytical: W=32, m=5', ...
    'Sim: W=128, m=3', 'Analytical: W=128, m=3' ...
});
grid on;
hold off;
