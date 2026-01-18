%% 2D Beamsteering function
%
% Inputs are in the shape of beamsteering(theta, num_ant, %optional% lambda)
% num_ant: số lượng anten trong mảng;  
% theta - góc tới (the angle of arrival) tính bằng radian
% lambda is a ratio of antenna spacing so với bước sóng
% 
% Default value is 1/2 corresponding to lambda/2 antenna spacing
%
% Returns beamsteering vector a(theta) expanded into the last dimension
% The output has the shape (size(theta) x num_ant)
% tạo vector beamsteering cho một mảng anten tuyến tính trong không gian 2D
function [a_theta] = beamsteering(varargin) 
%% - Nếu chỉ có 2 tham số → lấy theta và num_ant, còn lambda mặc định là 1/2.
%    Nếu có 3 tham số → lấy thêm lambda.
    if nargin == 2
        theta = varargin{1};
        num_ant = varargin{2};
    end
    if nargin > 2
        lambda = varargin{3};
    else
        lambda = 1/2;
    end
%% công thức
    dim = length(size(theta));                           %% só chiều của theta
    k = reshape([0:num_ant-1], [ones(1, dim), num_ant]); %% chỉ số anten
    
    a_theta = exp(1j*2*pi*lambda*k.*cos(theta));
%% Đầu ra: a_theta : ma trận phức, mỗi hàng là vector beamsteering ứng với một góc.    
end

