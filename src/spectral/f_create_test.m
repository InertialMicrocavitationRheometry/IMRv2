clc; clear all; close all;

addpath('../goldendata');
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
%find number of random characters to choose from
numRands = length(s); 
%specify length of random string to generate
sLength = 7;
%generate random string
file = s( ceil(rand(1,sLength)*numRands) );
filename = strcat('../goldendata/',file); 
[t,R,U] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',1,...
    'vapor',0,'plot',0);
display(filename)
save(filename,"t","R","U")

file = s( ceil(rand(1,sLength)*numRands) );
filename = strcat('../goldendata/',file); 
[t,R,U] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',2,...
    'vapor',0,'plot',0);
display(filename)
save(filename,"t","R","U")

file = s( ceil(rand(1,sLength)*numRands) );
filename = strcat('../goldendata/',file); 
[t,R,U] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',3,...
    'vapor',0,'plot',0);
display(filename)
save(filename,"t","R","U")

file = s( ceil(rand(1,sLength)*numRands) );
filename = strcat('../goldendata/',file); 
[t,R,U] = f_imrv2('r0',10E-6,'req',1E-6,'tvector',[0 5E-6],'radial',4,...
    'vapor',0,'plot',0);
display(filename)
save(filename,"t","R","U")