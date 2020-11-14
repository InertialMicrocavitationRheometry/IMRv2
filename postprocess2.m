%%  Original Post Process Deval Code

function [tsol_new, rsol_new, rdotsol_new] = postprocess2(tsol,rsol,rdotsol,sol)

    sign_change = find(diff(sign(rdotsol)));


    %rfill = [];
    %trange = linspace(0,tfinal,2000);
    %tfill = [];
    %trange = ODE.x;
    %idx_fill = 0;
    %t_new = zeros(1,length(trange)+Nfill*length(sign_change));
    Nfill = 3000;
    ishift = 0;
    t_new = tsol;

    for i = 1:length(sign_change)

    
    
        ttest = linspace(tsol(sign_change(i)), tsol(sign_change(i)+1),Nfill);
        %ytest = deval(sol,ttest);
        idx_beg = sign_change(i) + ishift;
        idx_end = idx_beg + 1;
    
    
        t_new = [t_new(1:idx_beg) ttest t_new(idx_end:end)];
    
        ishift = ishift + Nfill;
    
        %rfill = [rsol(1:idx_beg); ytest(1,:); rsol(1,idx_end,end)];
    
    
    end
    
    tsol_new = t_new;
    solve_new = deval(sol,tsol_new);
    rsol_new = solve_new(1,:);
    rdotsol_new = solve_new(2,:);


end 