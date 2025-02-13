clc; clear all; close all;

addpath('../src');
addpath('../unit_tests/');
load('file_ids.mat');

% equation options
tvector = linspace(0,20E-6,200);
count = 1;
for radial = 1:4
    for vapor = 0:1
        for bubtherm = 0:1
            for medtherm = 1
            filename = strcat('../unit_tests/',ids{count},'.dat'); 
            varin = {'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
                'vapor',vapor,'medtherm',medtherm};
            [ts,Rs,Us] = m_imrv2_spectral(varin{:},'Nt',15);
            [tf,Rf,Uf] = m_imrv2_finitediff(varin{:},'Nt',50);
            display(filename)
            if (norm(Rs-Rf,2) < 1E-2)
                disp('----> SUCCESS! <------')
                save(filename,"ts","Rs","Us")
            else
                disp('error radial not working')
                figure(count)
                hold on;
                plot(ts,Rs,'r--s')
                plot(tf,Rf,'k-.^')
            end
            count = count + 1;
            end
        end
    end
end