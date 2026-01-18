% Radar LoS channel generation
% The outer product of LoS beamsteering vectors are multipled with RCS 
% kênh radar trong điều kiện line of sight
% Trong radar, tín hiệu phát từ anten TX đi đến mục tiêu, phản xạ lại (Radar Cross Section – RCS), rồi được anten RX thu về.
%  xây dựng tensor kênh radar LoS bằng cách kết hợp:
%1- Beamsteering vector từ TX → Target.
%2- Beamsteering vector từ Target → RX.
%3- Radar Cross Section (RCS) của mục tiêu.

function [H] = radar_LOS_channel(target_pos, TX_pos, RX_pos, TX_ant, RX_ant, sigmasq_RCS)

    M_t = size(TX_pos, 1);
    M_r = size(RX_pos, 1);
    T = size(target_pos, 1);
% 1    TX đến target
    [angle, ~] = compute_angle_dist(TX_pos, target_pos);
    H = beamsteering(angle.', TX_ant);
    H1 = reshape(H, [T, M_t, 1, TX_ant, 1]);
% 2    target về TX
    [angle, ~] = compute_angle_dist(target_pos, RX_pos);
    H = beamsteering(angle, RX_ant);
    H2 = reshape(H, [T, 1, M_r, 1, RX_ant]);
% 3   sinh ngẫu nhiên - khả năng phản xạ của mục tiêu
    RCS = randn(T, M_t, M_r, 1, 1) .* sqrt(sigmasq_RCS);
    H = H1.*H2.*RCS;

end