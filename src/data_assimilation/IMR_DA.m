%

%
function [x,E2,Model_params_post,params]= IMR_DA( Model_params, load_info)
    
    % data import
    
    % 'sim' or 'exp'
    data_type      =   load_info{1};
    data_set       =   load_info{2};
    data_filepath  =   load_info{3};
    % name of file containing R vs T data
    data_filename  =   load_info{4};
    % number of 'peaks' to assimilate in radius data
    % (peaks = collapse points as in Estrada paper)
    num_peaks      =   2;
    
    
    import_data;
    
    % run main for corresponding method:
    [x,E2,P_post,params]  = main_En4D_peaks_par(Model_params,t,yth,R0_all,Req_all,tspan_all,peak_time_idx);
    
    %  collecting outputs
    Model_params_post{1}  =  Model_params{1};
    Model_params_post{2}  =  Model_params{2};
    Model_params_post{3}  =  P_post;
    
end
