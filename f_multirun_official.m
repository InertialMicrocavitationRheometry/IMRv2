function [] = f_multirun_official()

%% Clear everything 
close all; clear all; clc;
%%
N = 8^4;
addpath('./src/spectral/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One parameter cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newtonian model
mu_range1 = logspace(-6,-1,N);
parfor i = 1:N
    [t,R] = f_imrv2('voigt',1,'mu',mu_range1(i),'g',0)
    newtonian{i} = [t,R];
end
% Neo-Hookean model
G_range1 = logspace(-6,-1,N);
parfor i = 1:N
    disp(i/N)
    [t,R] = f_imrv2('neohook',1,'g',G_range1(i))
    neohook{i} = [t,R];
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two parameter case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kelvin-Voigt model
mu_range2 = logspace(-6,-1,N^(1/2));
G_range2 = logspace(-6,-1,N^(1/2));

% Create grid
[mu2, G2] = ndgrid(mu_range2, G_range2);
gridPoints2 = cell(length(mu_range2), length(G_range2));

% Fill in grid
for i = 1:length(mu_range2)
    for j = 1:length(G_range2)
        gridPoints2{i,j} = [mu2(i,j), G2(i,j)];
    end
end
cell2param = reshape(gridPoints2,[1,numel(gridPoints2)]);
parfor i = 1:N
    percent = i/N
    [t,R] = f_imrv2('voigt',1,'mu',cell2param{i}(1),'g',cell2param{i}(2));
    voigt{i} = [t,R];
end
voigt = reshape(voigt,[length(mu_range2), length(G_range2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three parameter case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upper-convected Maxwell model
mu_range3 = logspace(-6,-1,ceil(N^(1/3)));
G_range3 = logspace(-6,-1,ceil(N^(1/3)));
lambda1_range3 = logspace(-6,-1,ceil(N^(1/3)));
[mu3, G3, lambda13] = ndgrid(mu_range3,G_range3,lambda1_range3);
gridPoints3 = cell(length(mu_range3), length(G_range3), length(lambda1_range3));

for i = 1:length(mu_range3)
    for j = 1:length(G_range3)
        for k = 1:length(lambda1_range3)
            gridPoints3{i,j,k} = [mu3(i,j,k), G3(i,j,k), lambda13(i,j,k)];
        end
    end
end
cell3param = reshape(gridPoints3,[1,numel(gridPoints3)]);
parfor i = 1:N
    percent = i/N
    [t,R] = f_imrv2('oldb',1,'mu',cell3param{i}(1),'g',cell3param{i}(2),...
        'lambda1',cell3param{i}(3)*cell3param{i}(1)/cell3param{i}(2),'lambda2',0);
    ucm{i} = [t,R];
end
ucm = reshape(ucm,[length(mu_range3), length(G_range3), length(lambda1_range3)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Four parameter case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oldroyd-B model
mu_range4 = logspace(-6,-1,N^(1/4));
G_range4 = logspace(-6,-1,N^(1/4));
lambda1_range4 = logspace(-6,-1,N^(1/4));
lambda2_range4 = logspace(-6,-1,N^(1/4));
[mu4, G4, lambda14, lambda24] = ndgrid(mu_range4,G_range4,lambda1_range4, ...
                lambda2_range4);
gridPoints4 = cell(length(mu_range4), length(G_range4), length(lambda1_range4), ...
    length(lambda2_range4));

for i = 1:length(mu_range4)
    for j = 1:length(G_range4)
        for k = 1:length(lambda1_range4)
            for l = 1:length(lambda2_range4)
                gridPoints4{i,j,k,l} = [mu4(i,j,k,l), G4(i,j,k,l), ...
                    lambda14(i,j,k,l), lambda24(i,j,k,l)];
            end
        end
    end
end
cell4param = reshape(gridPoints4,[1,numel(gridPoints4)]);
parfor i = 1:N
    percent = i/N
    [t,R] = f_imrv2('oldb',1,'mu',cell4param{i}(1),'g',cell4param{i}(2),...
        'lambda1',cell4param{i}(3)*cell4param{i}(1)/cell4param{i}(2),...
        'lambda2',cell4param{i}(4)*cell4param{i}(1)/cell4param{i}(2));
    oldb{i} = [t,R];
end
oldb = reshape(oldb,[length(mu_range4), length(G_range4), ...
    length(lambda1_range4), length(lambda2_range4)]);



end


