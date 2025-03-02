clc;
clear;
close;

addpath('../toolchain/');
addpath('../src');
addpath('../tests/');
load('file_ids.mat');

% equation options
tvector = linspace(0,15E-6,100);
count = 1;
for radial = 1:4
    for vapor = 0:1
        for bubtherm = 0:1
            for medtherm = 0:1
                for stress = 0:5
                    filename1 = strcat('../tests/',ids{count+0},'.mat');
                    filename2 = strcat('../tests/',ids{count+1},'.mat');
                    varin = {'radial',radial,...
                        'bubtherm',bubtherm,...
                    'tvector',tvector,...
                        'vapor',vapor,...
                    'medtherm',medtherm,...
                        'stress',stress};
                    [tf,Rf] = m_imr_finitediff(varin{:},'Nt',100,'Mt',100);
                    [ts,Rs] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
                    if (norm(Rs-Rf,2) < 1E-2)
                        disp('----> SUCCESS! <------');
                        save(filename1,"Rf");
                        save(filename2,"Rs");
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

masstrans = 1;
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for vapor = 0:1
                for stress = 0:5
                    filename1 = strcat('../tests/',ids{count+0},'.mat');
                    varin = {'radial',radial,...
                        'bubtherm',bubtherm,...
                    'tvector',tvector,...
                        'vapor',vapor,...
                    'medtherm',medtherm,...
                        'stress',stress,...
                    'masstrans',masstrans};
                    [~,Rf] = m_imr_finitediff(varin{:},'Nt',100,'Mt',100);
                    save(filename1,"Rf");
                    count = count + 1;
                end
            end
        end
    end
end

masstrans = 0;
Req = 100e-6;
R0 = 100e-6;
T8 = 293;
tfin = 1E-3;
tvector = linspace(0,tfin,100);
for radial = 1:4
    for vapor = 0:1
        for bubtherm = 0:1
            for medtherm = 0:1
                for stress = 0:5
                    filename1 = strcat('../tests/',ids{count+0},'.mat');
                    filename2 = strcat('../tests/',ids{count+1},'.mat');
                    varin = {'radial',radial,...
                        'bubtherm',bubtherm,...
                    'tvector',tvector,...
                        'vapor',vapor,...
                    'medtherm',medtherm,...
                        'stress',stress,...
                    'Req',Req,...
                        'R0',R0,...
                    't8',T8};
                    [~,Rf] = m_imr_finitediff(varin{:},'Nt',100,'Mt',100);
                    [~,Rs] = m_imr_spectral(varin{:},'Nt',12,'Mt',12);
                    if (norm(Rf,2)/10-1 < 1e-6 && norm(Rs,2)/10-1 < 1e-6 )
                        disp('----> SUCCESS! <------');
                        save(filename1,"Rf");
                        save(filename2,"Rs");
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
