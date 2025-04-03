% file f_imr_da.m
% brief contains function f_imr_da for data assimilation

% brief This function conducts the IMR data assimilation routine to
% forecast material properties
function [x,ensemble,params_post,params]= m_imr_da(param_pre,load_info)
    
    addpath ../src;
    
    % data import
    data_filepath = load_info{3};
    % name of file containing R vs T data
    data_filename = load_info{4};
    % number of 'peaks' to assimilate in radius data
    % (peaks = collapse points as in Estrada paper)
    num_peaks = 2;
    
    % This file exports data in vector yth, it must be re-written for each data
    % set to ensure the data is formatted correctly
    
    % new data format
    exp_data0 = load([data_filepath,data_filename]);
    exp_data = exp_data0.sim_data;
    
    % number of time steps in experimental data
    n_exp = size(exp_data,2);
    
    peak_time_idx = zeros(size(n_exp));
    tspan_all = zeros(n_exp,1);
    t_exp_length = length(exp_data(1).t_norm(:));
    t_exp = zeros(n_exp,t_exp_length);
    R0_all = zeros(n_exp,1);
    max_index = 1;
    max_index_shift = 80;
    yth = zeros(n_exp,max_index_shift);
    Req_all = zeros(n_exp,1);
    t = zeros(n_exp,max_index_shift);
    
    for exp_idx = 1:n_exp
        
        R_exp = exp_data(exp_idx).Roft(1:1:end);
        t_exp(exp_idx,:) = exp_data(exp_idx).t_norm(1:1:end);
        R0_all(exp_idx) = exp_data(exp_idx).R0;
        
        % fixing shape if needed
        if length(R_exp(1,:)) == 1
            R_exp = R_exp';
        end
        if length(t_exp(1,:)) == 1
            t_exp = t_exp';
        end
        
        t_exp(exp_idx,:) = t_exp(exp_idx,:)- (t_exp(exp_idx,1)*2-t_exp(exp_idx,2));
        yth(exp_idx,:) = R_exp(max_index:max_index+max_index_shift)./R0_all(exp_idx);
        t(exp_idx,:)   = t_exp(exp_idx,max_index:max_index+max_index_shift);
        
        Req_all(exp_idx) = exp_data(exp_idx).Req;
        t(exp_idx,:) = t(exp_idx,:)- (t(exp_idx,1));
        
        % delete NaNs
        kk = 1;
        for jj = 1:size(yth,2)
            if isnan(yth(exp_idx,jj))
            else
                yth(exp_idx,kk) = yth(exp_idx,jj);
                t(exp_idx,kk) = t(exp_idx,jj);
                kk = kk + 1;
            end
        end
        yth(exp_idx,:) = yth(exp_idx,1:kk-1);
        t(exp_idx,:)   = t(exp_idx,1:kk-1);
        
        tspan_all(exp_idx) = t(exp_idx,end)-t(exp_idx,1);
        
        % find peak_time
        peak_indices = find(islocalmin(mean(yth(:,:),1)));
        peak_time_idx(exp_idx) = peak_indices(num_peaks);
    end
    
    % run main for corresponding method:
    [x,ensemble,p_post,params] = f_En4DVar(param_pre,t,yth,R0_all, ...
        Req_all,tspan_all,peak_time_idx);
    
    % collecting outputs
    params_post{1}  =  param_pre{1};
    params_post{2}  =  param_pre{2};
    params_post{3}  =  p_post;
    
end
