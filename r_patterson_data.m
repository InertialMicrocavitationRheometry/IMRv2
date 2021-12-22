
vals = readmatrix('./data/patterson_waveform_data.csv');

%Mostly Linear
t_ML = vals(:,2); %time [nondimensional]
p_ML = vals(:,3); %input pressure [nondimensional]

%Moderately Nonlinear
t_MN = vals(:,6); %time [nondimensional]
t_MN(isnan(t_MN)) = [];                
p_MN = vals(:,7); %input pressure [nondimensional]
p_MN(isnan(p_MN)) = [];
remove_index = [];
for i= 1:length(t_MN)-1
    if(t_MN(i) == t_MN(i+1) )
        remove_index = [remove_index i];
    end
end

for i = 1:length(remove_index)
    t_MN(remove_index(i)) = [];
    p_MN(remove_index(i)) = [];
end

%Highly Nonlinear
t_HN = vals(:,10); %time [nondimensional]
t_HN(isnan(t_HN)) = [];
p_HN = vals(:,11); %input pressure [nondimensional]
p_HN(isnan(p_HN)) = [];

remove_index = [];
for i= 1:length(t_HN)-1
    if(abs(t_HN(i) - t_HN(i+1)) < 0.01 )
    %if(t_HN(i) == t_HN(i+1))
        remove_index = [remove_index i];
    end
end

for i = 1:length(remove_index)
    t_HN(remove_index(i))
    t_HN(remove_index(i)) = [];
    p_HN(remove_index(i)) = [];
end
t_HN(89) = []; 
p_HN(89) = [];
remove_index = [remove_index 89]; %update remove index manually