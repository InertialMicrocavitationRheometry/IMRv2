clc;
clear;
close;

addpath('src');

% equation options
tvector = linspace(0,20E-6,200);
for radial = 1:4
    for bubtherm = 0:1
        for medtherm = 0:1
            for vapor = 0:1
                for stress = 0:5
                    for masstrans = 1
                        varin = {'radial',radial,...
                            'bubtherm',bubtherm,...
                            'tvector',tvector,...
                            'vapor',vapor,...
                            'medtherm',medtherm,...
                            'stress',stress,...
                            'masstrans',masstrans};
                        [tf,Rf,Uf] = m_imrv2_finitediff(varin{:},'Nt',20,'Mt',20);
                    end
                end
            end
        end
    end
end

% figure(1)
% hold on;
% plot(tf,Rf,'^')
