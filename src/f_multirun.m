function [vecoutput] = f_multirun()
N_Re = 100;
Re_vec = logspace(0,3,N_Re)
for i=0:1
    for j = 1:N_Re
        [t,R,p,U,trr,t00,I,T,C] = f_imrv2('Re', Re_vec(j))
    end
end
            
        