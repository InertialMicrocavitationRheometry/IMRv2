clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

% equation options
tvector = linspace(0,20E-6,200);
count = 1;
for radial = 1:4
    for vapor = 0:1
        for bubtherm = 0:1
            for medtherm = 0:1
                for stress = 0:4
                    filename1 = strcat('../tests/',ids{count+0},'.mat');
                    filename2 = strcat('../tests/',ids{count+1},'.mat');
                    varin = {'radial',radial,'bubtherm',bubtherm,'tvector',tvector,...
                        'vapor',vapor,'medtherm',medtherm,'stress',stress};
                    [tf,Rf,Uf] = m_imrv2_finitediff(varin{:},'Nt',100,'Mt',100);
                    [ts,Rs,Us] = m_imrv2_spectral(varin{:},'Nt',12,'Mt',12);
                    if (norm(Rs-Rf,2) < 1E-2)
                        disp('----> SUCCESS! <------');
                        save(filename1,"tf","Rf");
                        save(filename2,"ts","Rs");
                    else
                        disp('error radial not working')
                        figure(count)
                        hold on;
                        plot(ts,Rs,'r--s');
                        plot(tf,Rf,'k-.^');
                    end
                    count = count + 2;
                end
            end
        end
    end
end
