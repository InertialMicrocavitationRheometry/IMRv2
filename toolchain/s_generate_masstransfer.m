% file s_generate_masstransfer.m
% brief contains script to generate golden data for IMR mass transfer cases

% brief This script generates the golden data for cases using mass transfer
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

total_comb = 4*2*2*6*2*2*2*2*2;

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 1;
R0 = 2e-04;
Req = 3.5e-05;
radial_vec = 1:4;
bubtherm_vec = 0:1;
medtherm_vec = 0:1;
stress_vec = 0:5;
count = 1;

% precompute combinations
total_combinations = 4*2*2*6;

filenames_mass = cell(total_combinations,1);
for idx = 1:total_combinations
    filenames_mass{idx} = sprintf('../tests/%s.mat', ids{idx+total_comb});
end

% ensure the pool is active
if isempty(gcp('nocreate'))
    parpool('local', 8); % or whatever size you want
end

% set up combinations
dims = [length(stress_vec), length(medtherm_vec), ...
    length(bubtherm_vec), length(radial_vec)];
total_combinations = prod(dims);

% parallel dispatch
futures(total_combinations) = parallel.FevalFuture;

for idx_mass = 1:total_combinations
    [stress_idx, medtherm_idx, bubtherm_idx, radial_idx] = ind2sub(dims, idx_mass);
    
    radial   = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress   = stress_vec(stress_idx);
    
    varin = {'radial', radial, ...
        'bubtherm', bubtherm, ...
        'tvector', tvector, ...
        'vapor', vapor, ...
        'medtherm', medtherm, ...
        'stress', stress, ...
        'collapse', collapse, ...
        'r0', R0, ...
        'req', Req, ...
        'masstrans', masstrans};
    
    futures(idx_mass) = parfeval(@m_imr_fd_wrapper, 1, varin);
end

% save results in the correct order
for idx_mass = 1:total_combinations
    [completedIdx, result] = fetchNext(futures);
    savefile_fd(filenames_mass{completedIdx}, result);
    fprintf('Saved result %d to %s\n', completedIdx, filenames_mass{completedIdx});
end

delete(gcp('nocreate'));

% deterministic setup per worker
function Rm = m_imr_fd_wrapper(varin)
    [~, Rm] = m_imr_fd(varin{:}, 'Nt', 70, 'Mt', 70);
end

% savefile function
function savefile_fd(filename, data)
    Rm = data;
    save(filename, 'Rm');
end
