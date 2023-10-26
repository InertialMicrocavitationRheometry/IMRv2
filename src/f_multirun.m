function [vecoutput] = f_multirun()
close all; clear all; clc;
%consider neoHook and linear elastic 
%for Newtonian and elastic materials
%viscosity will range from zero to 10^1
%shear modulus will range from zero to 10^3

%starting with linear elasticity
%linear elasticity = 1 ; neoHook = 0
muvec = [0,logspace(-6,1,10)];
Gvec = [0,logspace(-3,3,10)];
for i=1:length(muvec)
    for j = 1:length(Gvec)
        mu = muvec(i);
        G = Gvec(j);
        filename = strcat('data_',num2str(mu),'_',num2str(G),'.mat');
        [t,R,p,Rdot,trr,t00,~,T,C,TL] = f_imrv2('linelas',1,'neoHook',0,...
         'mu',mu,'g',Gu);
    save(filename);
    end
end

%neoHook = 1; linear elasticity =0
for i=1:length(muvec)
    for j = 1:length(Gvec)
        mu = muvec(i);
        G = Gvec(j);
    [t,R,p,Rdot,trr,t00,~,T,C,TL] = f_imrv2('linelas',0,'neoHook',1,...
        'mu',mu,'g',Gu);
    save
    end
end
end 
        