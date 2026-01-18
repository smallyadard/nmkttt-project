%% tối đa hóa tín hiệu mong muốn tại đầu thu
% - beamforming vector được chọn trùng với kênh truyền H.
function F = beam_matched_filter(H)
    F = H;
    F = beam_normalization(F, 'UE');
end
     