
function [tsol_new, rsol_new, rdotsol_new] = postprocess(tsol,rsol,rdotsol,sol)

    %CODE TO FILL IN VECTOR
    sign_change = find(diff(sign(rdotsol)));
    Nfill = 3000;
    ishift = 0;
    a = 2;
    t_new = tsol;

    for i = 1:length(sign_change)
            
            if(rdotsol(sign_change(i)) > rdotsol(sign_change(i)+1))
                rg = rdotsol(sign_change(i));
                rl = rdotsol(sign_change(i)+1);
            else
                rg = rdotsol(sign_change(i)+1);
                rl = rdotsol(sign_change(i));
            end
        
            deltat = tsol(sign_change(i)+1)-tsol(sign_change(i));
            
            tstart = tsol(sign_change(i)+1) - deltat/((rg/rl)+1);
            tend = tsol(sign_change(i)+1);
            t_add = linspace(tstart,tend,Nfill);
            
            idx_beg = sign_change(i) + ishift;
            idx_end = idx_beg + 1;
        
            t_new = sort([t_new(1:idx_beg) t_add t_new(idx_end:end)]);
    
            ishift = ishift + Nfill;
    end
    remove_tnew = t_new > sol.x(end);
    t_new(remove_tnew) = [];

    tsol_new = t_new;
    solve_new = deval(sol,tsol_new);
    rsol_new = solve_new(1,:);
    rdotsol_new = solve_new(2,:);


end

%{  
FIRST VERSION OF CODE
%CODE TO FILL IN VECTOR
    sign_change = find(diff(sign(rdotsol)));
    Nfill = 3000;
    ishift = 0;
    a = 2;
    t_new = tsol;

    for i = 1:length(sign_change)

            %tm = ((tsol(sign_change(i)) - tsol(sign_change(i)+1)) / ...
            %(rdotsol(sign_change(i)) - rdotsol(sign_change(i)+1)))...
            %* -rdotsol(sign_change(i)+1) + tsol(sign_change(i)+1);
            
            tm = interp1(rdotsol((sign_change(i)-5):(sign_change(i)+5)),...
                    tsol((sign_change(i)-5):(sign_change(i)+5)),0);
            deltat = (tsol(sign_change(i))+tm)/a;
            t0 = tm - deltat;
            tf = tm + deltat;        
            t_add = linspace(t0,tf,Nfill);
        
        
            idx_beg = sign_change(i) + ishift;
            idx_end = idx_beg + 1;
        
            t_new = sort([t_new(1:idx_beg) t_add t_new(idx_end:end)]);
    
            ishift = ishift + Nfill;
    end
    remove_tnew = t_new > sol.x(end);
    t_new(remove_tnew) = [];

    tsol_new = t_new;
    solve_new = deval(sol,tsol_new);
    rsol_new = solve_new(1,:);
    rdotsol_new = solve_new(2,:);
%}

