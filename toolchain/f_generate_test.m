clc; clear all; close all;

addpath('../src');
addpath('../unit_tests/');
load('file_ids.mat');

% equation options
testb = 1;
teste = 4;
% radial tests
for id = testb:teste
    filename = strcat('../unit_tests/',ids{id},'.dat'); 
    [ts,Rs,Us] = m_imrv2_spectral('radial',id);
    [tf,Rf,Uf] = m_imrv2_finitediff('radial',id);
    display(filename)
    if (norm(Rs-Rf,2) < 1E-15)
        disp('----> SUCCESS! <------')
        save(filename,"ts","Rs","Us")
    else
        disp('error radial not working')
    end
end

% bubtherm tests
tvector = linspace(0,20E-6,200);
testb = testb+teste;
teste = teste + 4;
k = 1;
for id = testb:teste
    filename = strcat('../unit_tests/',ids{id},'.dat'); 
    [ts,Rs,Us] = m_imrv2_spectral('radial',k,'bubtherm',1,'Nt',15,'tvector',tvector);
    [tf,Rf,Uf] = m_imrv2_finitediff('radial',k,'bubtherm',1,'Nt',50,'tvector',tvector);
    display(filename)
    if (norm(Rs-Rf,2) < 1E-2)
        disp('----> SUCCESS! <------')
        save(filename,"ts","Rs","Us")
    else
        disp('error radial not working')
        figure(k)
        hold on;
        plot(ts,Rs,'r--s')
        plot(tf,Rf,'k-.^')
    end
    k = k + 1;
end

% bubtherm + vapor tests
tvector = linspace(0,20E-6,200);
testb = testb+teste;
teste = teste + 4;
k = 1;
for id = testb:teste
    filename = strcat('../unit_tests/',ids{id},'.dat'); 
    [ts,Rs,Us] = m_imrv2_spectral('radial',k,'bubtherm',1,'Nt',15,'tvector',tvector,'vapor',1);
    [tf,Rf,Uf] = m_imrv2_finitediff('radial',k,'bubtherm',1,'Nt',50,'tvector',tvector,'vapor',1);
    display(filename)
    if (norm(Rs-Rf,2) < 1E-2)
        disp('----> SUCCESS! <------')
        save(filename,"ts","Rs","Us")
    else
        disp('error radial not working')
        figure(k)
        hold on;
        plot(ts,Rs,'r--s')
        plot(tf,Rf,'k-.^')
    end
    k = k + 1;
end

%%
testb = testb + 1;
teste = teste + 1;
% medtherm tests
for id = testb:teste
    filename = strcat('../unit_tests/',ids{id},'.dat'); 
    [t,R,U] = m_imrv2_spectral('medtherm',1);
    display(filename)
    save(filename,"t","R","U")
end

% testb = testb + 5;
% teste = teste + 5;
% stressvc = 1:5;
% j = 0;
% % stress tests
% for id = testb:teste
%     j = j + 1;
%     filename = strcat('../unit_tests/',ids{id},'.dat'); 
%     [t,R,U] = m_imrv2_spectral('stress',stressvc(j));
%     display(filename)
%     save(filename,"t","R","U")
% end

% testb = testb + 1;
% teste = teste + 1;
% eps = [0.1,0.25];
% j = 0;
% % eps tests
% for id = testb:teste
%     j = j + 1;
%     filename = strcat('../unit_tests/',ids{id},'.dat'); 
%     [t,R,U] = m_imrv2_spectral('eps',epsvec(id));
%     display(filename)
%     save(filename,"t","R","U")
% end

testb = testb + 1;
teste = teste + 1;
% vapor tests
for id = testb:teste
    filename = strcat('../unit_tests/',ids{id},'.dat'); 
    [t,R,U] = m_imrv2_spectral('vapor',1);
    display(filename)
    save(filename,"t","R","U")
end
