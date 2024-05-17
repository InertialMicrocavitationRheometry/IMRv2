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
mu_range = logspace(mu_beg,mu_end,N);
G_range = logspace(G_beg,G_end,N);
lambda1_range = logspace(lambda1_beg,lambda1_end,N);
lambda2_range = logspace(lambda2_beg,lambda2_end,N);

% simulation count matrix to calculate the total number of simulations
sim_count_matrix = [1, 0, 0, 0;...
                    0, 1, 0, 0;...
                    1, 1, 0, 0;...
                    1, 1, 1, 0;...
                    1, 1, 1, 1];
% total number of simulations
num_sims = sum(N.^sum(sim_count_matrix,2));

% create the simulation parameter matrix
sim_param_matrix = zeros(num_sims,num_model_params);

sim_param_matrix = [
    table2array(combinations(mu_range)), zeros(N,3);    
    zeros(N,1), table2array(combinations(G_range)), zeros(N,2);    
    table2array(combinations(mu_range,G_range)), zeros(N^2,2);
    table2array(combinations(mu_range,G_range,lambda1_range)), zeros(N^3,1);
    table2array(combinations(mu_range,G_range,lambda1_range,lambda2_range))
    ];
