%% tính ma trận beamforming theo phương pháp Regularized Zero-Forcing (RZF) cho hệ thống MIMO
function F = beam_regularized_zeroforcing(varargin)
%%- H: mảng 3 chiều là kênh truyền có kích thước (U, M, N)
% P_comm: công suất truyền cho kênh thông tin.
% sigmasq_comm: phương sai nhiễu (noise variance).

    H = varargin{1};
    P_comm = varargin{2};
    sigmasq_comm = varargin{3};
    [U, M, N] = size(H);

    if nargin == 3
        lambda = U*sigmasq_comm/P_comm;
    else
        lambda = varargin{4};
    end
% công thức theo phương pháp Regularized Zero-Forcing (RZF)
    H_stacked = reshape(H, U, []);
    F_stacked = pinv(eye(U)*lambda + H_stacked * H_stacked')*H_stacked;
    F = reshape(F_stacked, U, M, N);
    F = beam_normalization(F, 'UE');
%% out: ma trận beamforming đã được chuẩn hóa, có cùng kích thước vói H .
end
     