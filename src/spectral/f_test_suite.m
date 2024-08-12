clc; clear all; close all;

addpath('../goldendata');
%find number of random characters to choose from
filename = strcat('../goldendata/','0VwpAlX'); 
load(filename);
[~,Rtest,Utest] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',1,...
    'vapor',0,'plot',0);
deltaR = sum(abs(R-Rtest));
deltaU = sum(abs(U-Utest));
if deltaR > 1E-12
    error('TEST 1 FAILS RAD: not within tolerance');
elseif deltaU > 1E-12
    error('TEST 1 FAILS VEL: not within tolerance');
else
    disp('TEST 1 PASSED');
end

filename = strcat('../goldendata/','4AcacvT'); 
load(filename);
[~,Rtest,Utest] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',2,...
    'vapor',0,'plot',0);
deltaR = sum(abs(R-Rtest));
deltaU = sum(abs(U-Utest));
if deltaR > 1E-12
    error('TEST 2 FAILS RAD: not within tolerance');
elseif deltaU > 1E-12
    error('TEST 2 FAILS VEL: not within tolerance');
else
    disp('TEST 2 PASSED');
end

filename = strcat('../goldendata/','wdCKsdJ'); 
load(filename);
[~,Rtest,Utest] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',3,...
    'vapor',0,'plot',0);
deltaR = sum(abs(R-Rtest));
deltaU = sum(abs(U-Utest));
if deltaR > 1E-12
    error('TEST 3 FAILS RAD: not within tolerance');
elseif deltaU > 1E-12
    error('TEST 3 FAILS VEL: not within tolerance');
else
    disp('TEST 3 PASSED');
end

filename = strcat('../goldendata/','VlLtP4Q'); 
load(filename);
[~,Rtest,Utest] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',4,...
    'vapor',0,'plot',0);
deltaR = sum(abs(R-Rtest));
deltaU = sum(abs(U-Utest));
if deltaR > 1E-12
    error('TEST 4 FAILS RAD: not within tolerance');
elseif deltaU > 1E-12
    error('TEST 4 FAILS VEL: not within tolerance');
else
    disp('TEST 4 PASSED');
end

clear;