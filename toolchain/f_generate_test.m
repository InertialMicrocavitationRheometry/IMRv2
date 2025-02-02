clc; clear all; close all;
addpath('src');
addpath('src/spectral/');
addpath('toolchain/');
addpath('unit_tests/');
load('file_ids.mat');

% equation options

testb = 1;
teste = 4;
% radial tests
for id = testb:teste
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('radial',id);
    display(filename)
    save(filename,"t","R","U")
end

testb = testb + 1;
teste = teste + 1;
% bubtherm tests
for id = testb:teste
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('bubtherm',1);
    display(filename)
    save(filename,"t","R","U")
end

testb = testb + 1;
teste = teste + 1;
% medtherm tests
for id = testb:teste
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('medtherm',1);
    display(filename)
    save(filename,"t","R","U")
end

testb = testb + 5;
teste = teste + 5;
stressvc = 1:5;
j = 0;
% stress tests
for id = testb:teste
    j = j + 1;
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('stress',stressvc(j));
    display(filename)
    save(filename,"t","R","U")
end

testb = testb + 1;
teste = teste + 1;
eps = [0.1,0.25];
j = 0;
% eps tests
for id = testb:teste
    j = j + 1;
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('eps',epsvec(id));
    display(filename)
    save(filename,"t","R","U")
end

testb = testb + 1;
teste = teste + 1;
% vapor tests
for id = testb:teste
    filename = strcat('unit_tests/',ids{id},'.dat'); 
    [t,R,U] = f_imrv2('vapor',1);
    display(filename)
    save(filename,"t","R","U")
end
