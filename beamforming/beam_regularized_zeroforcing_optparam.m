%% chọn lambda cho RZF để tìm tham số đem lại hiệu suất tốt nhất
function F = beam_regularized_zeroforcing_optparam(H_comm, P_comm, F_sensing, P_sensing, sigmasq_ue, lambda_min, lambda_max, num_lambda)
% F_sensing ma trận beamforming cho sensing 
% P_sensing công suất 
    
% tính SINR
eval_fun = @(x) compute_SINR(H_comm, x, P_comm, F_sensing, P_sensing, sigmasq_ue);

    f_rzf_fun = @(x) beam_regularized_zeroforcing(H_comm, P_comm, sigmasq_ue, x);
% chọn lambda bằng tìm kiếm tuyến tính    
    opt_lambda = linear_parameter_search(eval_fun, f_rzf_fun, lambda_min, lambda_max, num_lambda);
    
    F = f_rzf_fun(opt_lambda);
    %% đầu ra ma trận beamforming F
end
     