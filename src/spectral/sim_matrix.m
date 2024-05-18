clear all; clc; close all;

% this is the simulation matrix for each material in the physical world,
% e.g., gelatin, agarose, etc.

N = 10; % number of values per model parameter
num_model_params = 4; 
num_model_param_range = N*ones(num_model_params,1); 

% manually change ranges as needed
mu_beg = -6;
mu_end = -1;
G_beg = -6;
G_end = -1;
lambda1_beg = -6;
lambda1_end = -1;
lambda2_beg = -6;
lambda2_end = -1;

% creating the logspaces for each of the material parameters
N = 8;
N1 = N^4; N2 = N^2; N3 = ceil(N^(4/3)); N4 = N;
mu_range1 = logspace(mu_beg,mu_end,N1);
mu_range2 = logspace(mu_beg,mu_end,N2);
mu_range3 = logspace(mu_beg,mu_end,N3);
mu_range4 = logspace(mu_beg,mu_end,N4);
G_range1 = logspace(G_beg,G_end,N1);
G_range2 = logspace(G_beg,G_end,N2);
G_range3 = logspace(G_beg,G_end,N3);
G_range4 = logspace(G_beg,G_end,N4);
lambda1_range3 = logspace(lambda1_beg,lambda1_end,N3);
lambda1_range4= logspace(lambda1_beg,lambda1_end,N4);
lambda2_range4 = logspace(lambda2_beg,lambda2_end,N4);


% simulation count matrix to calculate the total number of simulations
sim_count_matrix = [1, 0, 0, 0;...
                    0, 1, 0, 0;...
                    1, 1, 0, 0;...
                    1, 1, 1, 0;...
                    1, 1, 1, 1];
% total number of simulations
num_sims = 5*N1;

% create the simulation parameter matrix
sim_param_matrix = zeros(num_sims,num_model_params);

sim_param_matrix = [
    table2array(combinations(mu_range1)), zeros(N1,3);    
    zeros(N1,1), table2array(combinations(G_range1)), zeros(N1,2);    
    table2array(combinations(mu_range2,G_range2)), zeros(N1,2);
    table2array(combinations(mu_range3,G_range3,lambda1_range3)), zeros(N1,1);
    table2array(combinations(mu_range4,G_range4,lambda1_range4,lambda2_range4))
    ];

model4grid_param_unsort = sim_param_matrix(4096*3+1:4096*4,:);

model4grid_param_sort = zeros(3,N3,N3,N3);
for l = 1:3
    for i = 1:N3
        for j = 1:N3
            for k = 1:N3
                idx = k + (j-1)*N3 + (i-1)*N3*N3;
                model4grid_param_sort(l,i,j,k) = model4grid_param_unsort(idx,l);
            end
        end
    end
end