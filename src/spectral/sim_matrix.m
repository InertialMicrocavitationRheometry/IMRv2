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
mu_range1 = logspace(mu_beg,mu_end,10000);
mu_range2 = logspace(mu_beg,mu_end,1000);
mu_range3 = logspace(mu_beg,mu_end,100);
mu_range4 = logspace(mu_beg,mu_end,10);
G_range1 = logspace(G_beg,G_end,10000);
G_range2 = logspace(G_beg,G_end,1000);
G_range3 = logspace(G_beg,G_end,100);
G_range4 = logspace(G_beg,G_end,10);
lambda1_range3 = logspace(lambda1_beg,lambda1_end,100);
lambda1_range4= logspace(lambda1_beg,lambda1_end,10);
lambda2_range4 = logspace(lambda2_beg,lambda2_end,10);
N1 = 10000; N2 = 1000; N3 = 100; N4 = 10;

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
    table2array(combinations(mu_range1)), zeros(N1,3);    
    zeros(N1,1), table2array(combinations(G_range1)), zeros(N1,2);    
    table2array(combinations(mu_range2,G_range2)), zeros(N2^2,2);
    table2array(combinations(mu_range3,G_range3,lambda1_range3)), zeros(N3^3,1);
    table2array(combinations(mu_range4,G_range4,lambda1_range4,lambda2_range4))
    ];
