function [F_star, feasible, SSNR_opt] = opt_jsc_FA_hybrid( ...
        H_comm, sigmasq_comm, feasible_SINR, ...
        sensing_beamsteering, sensing_streams, ...
        sigmasq_radar_rcs, P_all)
% HYBRID FIREFLY ALGORITHM 
% Enhancement: Levy flight + Elite guidance
% Objective: maximize SSNR under SINR & Power constraints

%% ================= 1. Problem setup =================
[U, M, N] = size(H_comm);
MN = M * N;
num_streams = U + sensing_streams;
dim = num_streams * 2 * MN;

Lb = -sqrt(P_all) * ones(1, dim);
Ub =  sqrt(P_all) * ones(1, dim);

%% ================= 2. Hybrid FA parameters =================
nPop = 30;
MaxIter = 100;

beta0 = 1.0;
gamma = 1.0;

alpha0 = 0.6;
alpha_min = 0.05;

p_levy = 0.3;        % xác suất Lévy flight
elite_ratio = 0.1;  % % cá thể ưu tú

%% ================= 3. Initialization =================
X = Lb + (Ub - Lb) .* rand(nPop, dim);
Cost = zeros(nPop,1);

%% ================= 4. Main loop =================
for iter = 1:MaxIter

    alpha = alpha0 * (1 - iter/MaxIter) + alpha_min;

    %% ---- Evaluate population ----
    for i = 1:nPop
        f_vecs = decode_chromosome(X(i,:), MN, num_streams);

        F_comm = reshape(f_vecs(:,1:U), [U, M, N]);
        F_sensing = reshape(f_vecs(:,U+1:end), [sensing_streams, M, N]);

        metric = compute_metrics( ...
            H_comm, F_comm, sigmasq_comm, ...
            sensing_beamsteering, F_sensing, sigmasq_radar_rcs);

        violation_sinr  = max(0, feasible_SINR - metric.min_SINR);
        violation_power = max(0, sum(metric.power) - P_all);

        Cost(i) = -metric.SSNR ...
                  + 1e6 * violation_sinr ...
                  + 1e6 * violation_power;
    end

    %% ---- Sort fireflies ----
    [Cost, idx] = sort(Cost);
    X = X(idx,:);

    best = X(1,:);

    %% ---- Move fireflies (Hybrid FA) ----
    for i = 2:nPop
        for j = 1:i-1
            r = norm(X(i,:) - X(j,:));
            beta = beta0 * exp(-gamma * r^2);

            % FA attraction
            X(i,:) = X(i,:) ...
                   + beta * (X(j,:) - X(i,:));
        end

        % Elite guidance
        X(i,:) = X(i,:) + 0.2 * (best - X(i,:));

        % Lévy flight (exploration)
        if rand < p_levy
            X(i,:) = X(i,:) + 0.01 * levy_flight(dim);
        end

        % Adaptive random walk
        X(i,:) = X(i,:) + alpha * (rand(1,dim) - 0.5);

        % Boundary control
        X(i,:) = max(X(i,:), Lb);
        X(i,:) = min(X(i,:), Ub);
    end
end

%% ================= 5. Output =================
best_x = X(1,:);
f_vecs = decode_chromosome(best_x, MN, num_streams);

F_star = zeros(MN, MN, num_streams);
for s = 1:num_streams
    F_star(:,:,s) = f_vecs(:,s) * f_vecs(:,s)';
end

F_comm = reshape(f_vecs(:,1:U), [U, M, N]);
F_sensing = reshape(f_vecs(:,U+1:end), [sensing_streams, M, N]);

metric = compute_metrics( ...
    H_comm, F_comm, sigmasq_comm, ...
    sensing_beamsteering, F_sensing, sigmasq_radar_rcs);

SSNR_opt = metric.SSNR;
feasible = (metric.min_SINR >= feasible_SINR) ...
           && (sum(metric.power) <= P_all);

end

%% ================= Auxiliary functions =================
function f_vecs = decode_chromosome(x, MN, num_streams)
    f_vecs = zeros(MN, num_streams);
    idx = 1;
    for s = 1:num_streams
        re = x(idx:idx+MN-1); idx = idx + MN;
        im = x(idx:idx+MN-1); idx = idx + MN;
        f_vecs(:,s) = re.' + 1i*im.';
    end
end

function step = levy_flight(d)
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
            (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    step = u ./ abs(v).^(1/beta);
end
