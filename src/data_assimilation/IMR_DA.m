function [x,E2,model_params_post,params]= IMR_DA( model_params )
    
    % data import
    
    % data_type      =   load_info{1}; % 'sim' or 'exp'
    % data_set       =   load_info{2};
    % data_filepath  =   load_info{3};
    % data_filename  =   load_info{4}; % name of file containing R vs T data
    % num_peaks      =   2; % number of 'peaks' to assimilate in radius data
    % (peaks = collapse points as in Estrada paper)
    import_data;
    
    % run main for corresponding method:
    
    [x,E2,P_post,params]  = main_En4D_peaks_par(model_params,t,yth,R0_all,Req_all,tspan_all,peak_time_idx);
    
    % collecting outputs
    
    model_params_post{1}  =  model_params{1};
    model_params_post{2}  =  model_params{2};
    model_params_post{3}  =  P_post;
    
end
