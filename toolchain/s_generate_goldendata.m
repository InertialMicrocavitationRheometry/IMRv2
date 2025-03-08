% file s_generate_goldendata.m
% brief contains script to generate golden data for IMR

% brief This script generates the golden data for various test suites used
% for IMR
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

% thermal test cases
tvector = linspace(0,15E-6,100);
threshold = 5e-2;
collapse = 1;
masstrans = 0;
vapor = 1;
count = 1;
counter = 1;
R0 = 50e-6;
Req = R0/12;
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for stress = 0:5
                filename1 = strcat('../tests/',ids{count+0},'.mat');
                filename2 = strcat('../tests/',ids{count+1},'.mat');
                varin = {'progdisplay',0,...
                    'radial',radial,...
                    'bubtherm',bubtherm,...
                    'tvector',tvector,...
                    'vapor',vapor,...
                    'medtherm',medtherm,...
                    'masstrans',masstrans,...
                    'collapse',collapse,...
                    'lambda2',0,...
                    'Req',Req,...
                    'R0',R0,...
                    'stress',stress};
                [tf,Rf] = m_imr_fd(varin{:},'Nt',50,'Mt',50);
                [ts,Rs] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
                if (norm(Rs-Rf,2) < threshold)
                    disp('----> SUCCESS! <------');
                    save(filename1,"Rf");
                    save(filename2,"Rs");
                else
                    figure(count)
                    hold on;
                    plot(ts,Rs,'r--s');
                    plot(tf,Rf,'k-.^');
                    error('error radial not working')
                end
                count = count + 2;
                counter = counter + 1;
            end
        end
    end
end

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 1;
R0 = 2e-04;
Req = 3.5e-05;
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for stress = 0:5
                filename1 = strcat('../tests/',ids{count+0},'.mat');
                varin = {'radial',radial,...
                    'bubtherm',bubtherm,...
                    'tvector',tvector,...
                    'vapor',vapor,...
                    'medtherm',medtherm,...
                    'stress',stress,...
                    'collapse',collapse,...
                    'r0',R0,...
                    'req',Req,...
                    'masstrans',masstrans};
                [~,Rf] = m_imr_fd(varin{:},'Nt',50,'Mt',50);
                save(filename1,"Rf");
                count = count + 1;
            end
        end
    end
end
