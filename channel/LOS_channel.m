% LOS channel from x_pos to y_pos with x_ant and y_ant=1
% sinh kênh truyền line of sight giữa điểm truy cập phat và thu, với nhiều anten phát và 1 anten thu
function [H] = LOS_channel(x_pos, y_pos, x_ant)
% tính góc tới giữa anten phát và anten thu
    [angle, ~] = compute_angle_dist(x_pos, y_pos); 
% tạo vector beamsteering tương ứng với số anten phát
    H = beamsteering(angle.', x_ant); % y by x
   
end