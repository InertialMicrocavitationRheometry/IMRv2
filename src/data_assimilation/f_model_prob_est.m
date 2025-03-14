% Evaluating model probability
% Inputs:  Q_model   --   All the available models with their posterior distributions
%          q_data    --   Measurements
function [model_prob_all] = f_model_prob_est(Q_model, q_data)
    [N,Nx]            =   size(q_data);
    N_model           =   size(Q_model,1);
    model_prob_all    =   zeros(N_model,1);
    % removed the constant factor here to avoid numerical issues
    norm_fun_log      =   @(mu, Sigma, x) -1/2*(+(((x-mu)/(Sigma))*(x-mu)'));
    
    % measurement variaces in the data
    sigma_data        =   zeros(Nx);
    mu_data           =   zeros(1,Nx);
    
    for j = 1:Nx
        [sigma_data(j,j),mu_data(j)] = robustcov(q_data(:,j));
    end
    
    P_EVD             =   zeros(N,1);
    for j1 = 1:N
        % data evidence
        P_EVD(j1)     =  (norm_fun_log(q_data(j1,:), sigma_data, mu_data))';
    end
    
    % calculate the marginal likelihood
    for k = 1:N_model
        q_model              =  Q_model{k};
        [M,~]                =  size(q_model);
        P_LKH                =  zeros(1,M);
        for j2 = 1:M
            for j1 = 1:N
                % marginal likelihood
                P_LKH(1,j2)  =  P_LKH(1,j2) + 1/N*(norm_fun_log(q_data(j1,:), ...
                    sigma_data, q_model(j2,:) ))';
            end
        end
        model_prob_all(k)    =  log10(mean(exp(P_LKH)));
    end
    % normalization
    model_prob_all         =  (10.^model_prob_all)./sum(10.^model_prob_all);
end
