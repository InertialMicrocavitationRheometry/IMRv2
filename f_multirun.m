function [vecoutput] = f_multirun()
close all; clear all; clc;
%consider neoHook and linear elastic 
%for Newtonian and elastic materials
%viscosity will range from zero to 10^1
%shear modulus will range from zero to 10^3

addpath('spectral');

%starting with linear elasticity
%linear elasticity = 1 ; neoHook = 0
muvec = [0,logspace(-6,1,10)];
Gvec = [0,logspace(-3,3,10)];

% for i=1:length(muvec)
%     for j = 1:length(Gvec)
%         mu = muvec(i);
%         G = Gvec(j);
%         foldername = '../bdata/linearelasticity';
%         mkdir(foldername);
%         filename = strcat(foldername,'/data_',num2str(mu),'_',num2str(G),'.mat');
%         [t,R,p,Rdot,trr,t00,~,T,C,TL] = f_imrv2('linelas',1,'neoHook',0,...
%          'mu',mu,'g',G);
%         save(filename);
%     end
% end
%neoHook = 1; linear elasticity =0
for i=1:length(muvec)
    for j = 1:length(Gvec)
        mu = muvec(i);
        G = Gvec(j);
        foldername = '../bdata/neoHookean';
        mkdir(foldername);
        filename = strcat(foldername,'/data_',num2str(mu),'_',num2str(G),'.csv');
        [t,R,p,Rdot,trr,t00,~,T,C,TL] = f_imrv2('linelas',0,'neoHook',1,...
        'mu',mu,'g',G);
        a = [t,R];
        writematrix(a, filename);
    end
end
%%
forlegend = [];
addpath("../IMR_v2_matlab")
filenames = dir("../PolyAcry*.mat");
%initialize the variables inside for loop?
% t = zeros(size(filenames))';
% R1 = zeros(size(filenames))';
% R2 = zeros(size(filenames))';
%extract bubble radius and time for each file
for k = 1:size(filenames)
    filename = filenames(k).name;
    %alldata = importdata(filename);
    alldata = load(filename);
    data = alldata.BubblePlotSet;
    t = data.time';
    R1 = data.radius(1,:)';
    R2 = data.radius(2,:)';
    a = [t,R1,R2];
    %save necessary data as csv
    filename = filename(1:end-4);
    csvwrite(filename,a)
    %writematrix(a,filename)
    %writematrix(a, [filename '.csv'])
    
    %for legend
    forlegend{k} = [filename];
    
    
    %if value is NaN, remove?
    %should also be saving dimensionless numbers: ReB, De, CA
    %problems: t, R1, R2 should be saving data for EACH filename, but is getting overwritten
    
    %Plot time and radius for each camera
    figure(1)
    hold on
    plot(t,R1,'.')
    legend(forlegend,'Location','best')
    
    figure(2)
    hold on
    plot(t,R2,'.')
    legend(forlegend,'Location','best')
end
rmpath("../IMR_v2_matlab")
rmpath('spectral');
end 

