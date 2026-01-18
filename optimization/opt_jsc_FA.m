function [F_star, feasible, SSNR_opt] = opt_jsc_FA(H_comm, sigmasq_comm, feasible_SINR, sensing_beamsteering, sensing_streams, sigmasq_radar_rcs, P_all)
% OPT_JSC_FA: Tối ưu Beamforming ISAC bằng Firefly Algorithm (FA)
% Đánh giá cost dựa trên compute_metrics (power, SINR, SSNR)

    %% 1. Thiết lập thông số bài toán
    [U, M, N] = size(H_comm); % U users, M APs, N antennas
    MN = M * N;
    num_streams = U + sensing_streams; % Tổng số vector cần tìm (User + Sensing)

    %% 2. Thông số FA
    n = 30;                 % Số lượng đom đóm
    MaxGeneration = 100;    % Số vòng lặp
    alpha = 0.5;            % Randomness
    beta0 = 1.0;            % Attractiveness
    gamma_fa = 1.0;         % Absorption

    dim = num_streams * 2 * MN; % Kích thước biến tối ưu

    % Giới hạn tìm kiếm
    Lb = -sqrt(P_all) * ones(1, dim);
    Ub =  sqrt(P_all) * ones(1, dim);

    % Khởi tạo quần thể
    fireflies = Lb + (Ub - Lb) .* rand(n, dim);
    LightIntensity = zeros(n, 1);

    %% 3. Vòng lặp FA
    for iter = 1:MaxGeneration
        % Đánh giá cost
        for i = 1:n
            f_vecs = decode_chromosome(fireflies(i,:), MN, num_streams);

            % Tách thành F_comm và F_sensing
            F_comm = reshape(f_vecs(:,1:U), [U, M, N]);
            F_sensing = reshape(f_vecs(:,U+1:end), [sensing_streams, M, N]);

            % Tính metric
            metric_struct = compute_metrics(H_comm, F_comm, sigmasq_comm, ...
                                            sensing_beamsteering, F_sensing, sigmasq_radar_rcs);

            % Phạt nếu vi phạm SINR hoặc Power
            violation_sinr  = max(0, feasible_SINR - metric_struct.min_SINR);
            violation_power = max(0, sum(metric_struct.power) - P_all);
            % Cost = -SSNR + penalty
            LightIntensity(i) = -metric_struct.SSNR ...
                                + 1e6*violation_sinr ...
                                + 1e6*violation_power;
        end

        % Sắp xếp theo cost
        [LightIntensity, sort_idx] = sort(LightIntensity);
        fireflies = fireflies(sort_idx, :);

        % Di chuyển đom đóm
        for i = 1:n
            for j = 1:n
                if LightIntensity(j) < LightIntensity(i)
                    r = norm(fireflies(i,:) - fireflies(j,:));
                    beta = beta0 * exp(-gamma_fa * r^2);

                    fireflies(i,:) = fireflies(i,:) ...
                        + beta * (fireflies(j,:) - fireflies(i,:)) ...
                        + alpha * (rand(1, dim) - 0.5);

                    % Kẹp vào biên
                    fireflies(i,:) = max(fireflies(i,:), Lb);
                    fireflies(i,:) = min(fireflies(i,:), Ub);
                end
            end
        end

        % Giảm alpha
        alpha = alpha * 0.98;
    end

    %% 4. Kết quả đầu ra
    best_firefly = fireflies(1,:);
    f_vecs = decode_chromosome(best_firefly, MN, num_streams);

    F_star = zeros(MN, MN, num_streams);
    for s = 1:num_streams
        F_star(:,:,s) = f_vecs(:,s) * f_vecs(:,s)';
    end

    % Tính metric cuối cùng
    F_comm = reshape(f_vecs(:,1:U), [U, M, N]);
    F_sensing = reshape(f_vecs(:,U+1:end), [sensing_streams, M, N]);
    metric_struct = compute_metrics(H_comm, F_comm, sigmasq_comm, ...
                                    sensing_beamsteering, F_sensing, sigmasq_radar_rcs);

    SSNR_opt = metric_struct.SSNR;
    feasible = (metric_struct.min_SINR >= feasible_SINR) && (sum(metric_struct.power) <= P_all);

end

%% --- Hàm phụ ---
function f_vecs = decode_chromosome(x, MN, num_streams)
    f_vecs = zeros(MN, num_streams);
    idx = 1;
    for s = 1:num_streams
        re_part = x(idx : idx+MN-1); idx = idx+MN;
        im_part = x(idx : idx+MN-1); idx = idx+MN;
        f_vecs(:,s) = re_part.' + 1i*im_part.';
    end
end