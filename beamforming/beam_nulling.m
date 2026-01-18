%% Function to project sensing beam to null-space of communication Channels
% điều chỉnh vector beamforming của sensing sao cho nó nằm trong không gian trực giao với kênh truyền thông .
% mục tiêu loại bỏ nhiễu
function [F_sensing] = beam_nulling(H_comm, sensing_beamsteering)

    [T, M, N] = size(sensing_beamsteering);
    F_sensing = zeros(T, M, N);

    % Compute the nullspace projection at each AP
    % Tại mỗi anten phát, chiếu beam sensing vào null-space của kênh truyền thông để loại bỏ nhiễu
    for m=1:M
        H_null = squeeze(H_comm(:, m, :)).'; % NxU channel at each AP
        for t = 1:T
            F_sensing(t, m, :) = (eye(N) - H_null * pinv(H_null'* H_null) * H_null') * squeeze(sensing_beamsteering(t, m, :));
        end
    end

    % Normalize the power
    F_sensing = beam_normalization(F_sensing, 'UE');
end

