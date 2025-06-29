% file s_resolution_check.m
% brief contains script to check resolution capability of IMR
clc;
clear;
close;
format long;
    
    addpath('../../src/forward_solver/');
    
    % material parameter test cases
    tvector = linspace(0,100E-6,500);
    collapse = 1;
    masstrans = 1;
    vapor = 1;
    bubtherm = 1;
    medtherm = 1;
    R0 = 200e-6;
    Req = R0/8;
    mu = 1E-1;
    G = 1E3;
    alphax = 0.5;
    lambda1 = 1E-3;
    lambda2 = 0;
    radial = 3;
    stress = 3;
    
    % define parameter vector and enforce Mt/Nt = 4
    Nt_vec = 2.^(5:9);
    nNt = numel(Nt_vec);
    
    tvec = cell(nNt,1);
    Rvec = cell(nNt,1);
    
    % parallel pool
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', 5);
    end
    
    % futures to track jobs
    futures = parallel.FevalFuture.empty(nNt,0);
    
    % launch all jobs
    for iNt = 1:nNt
        Nt = Nt_vec(iNt);
        Mt = 4*Nt;  % enforce ratio
        varin = {'progdisplay',0,...
            'radial',radial,...
            'bubtherm',bubtherm,...
            'tvector',tvector,...
            'vapor',vapor,...
            'medtherm',medtherm,...
            'masstrans',masstrans,...
            'collapse',collapse,...
            'lambda2',0,...
            'Req',Req,...
            'R0',R0,...
            'mu',mu,...
            'G',G,...
            'alphax',alphax,...
            'lambda1',lambda1,...
            'stress',stress,...
            'Nt',Nt,...
            'Mt',Mt};
        
        futures(iNt) = parfeval(pool, @f_imr_fd, 2, varin{:});
    end
    
    % collect results
    for idx = 1:nNt
        [completedIdx, tf, Rf] = fetchNext(futures);
        tvec{completedIdx} = tf;
        Rvec{completedIdx} = Rf;
        fprintf('Completed %d of %d (Nt=%d, Mt=%d)\n', ...
            idx, nNt, Nt_vec(completedIdx), 4*Nt_vec(completedIdx));
    end
    
    save("fdres_Nt_ratio4.mat","tvec","Rvec");
    delete(gcp('nocreate'));
