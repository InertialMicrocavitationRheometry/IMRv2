clear all; clc; close all;
Nv = 200;
Lvec = [3,4,5,6,7,10,20,0.1,0.01,0.3,0.5];

for i = 1:length(Lvec)
    Lv = Lvec(i);
    preStressInt(Lv,Nv);
end

function cdd = preStressInt(L,N)
    
    % integral precomputations
    Lstr = ['L' num2str(L,2)];
    disp(Lstr);
    Lstr = strrep(Lstr,'.','p');
    if exist('StressIntStore.mat','file') ~= 0
        load('StressIntStore.mat','store');
        if isfield(store,Lstr) == 1
            if size(store.(Lstr),2) >= N, Nstart = 0;
                else
                Nstart = size(store.(Lstr),2) + 1;
                disp('Past integral precomputation not found in full, catching up ...');
            end
            else
            Nstart = 1;
            disp('Past integral precomputation not found, starting anew ...');
        end
        else
        Nstart = 1;
        store = struct;
        disp('Past integral precomputation not found, starting anew ...');
    end
    if Nstart ~= 0 % begin extended precomputation
        
        store.(Lstr)(Nstart:N) = StressInt(L,N,Nstart);
        save('StressIntStore.mat','store');
        disp('Precomputation completed.');
        
    end
    cdd = store.(Lstr)(1:N)';
    
end


function cdd = StressInt(L,N,varargin)
    
    if nargin == 2
        k = 1;
        else
        k = varargin{1};
    end
    
    syms x;
    cdd = zeros(N-k+1,1);
    
    for n = k:N
        cdd(n-k+1) = subs(2*L*int((cos(n*acos(x))-1)/((L*(2/(1-x)-1)+1)*(1-x)^2),-1,1));
    end
    
end
