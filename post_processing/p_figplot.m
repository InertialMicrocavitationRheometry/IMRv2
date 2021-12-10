figure(1)
hold on;
% tspan = Tout*fnatural;
% Roft = (Rout(:,1)/Ro_w);
plot(t,R, lm,'LineWidth',2); 

% figure(2)
% gamma = -4*(Rout(:,2))./(Rout(:,1));
% gamman = gamma/fnatural;
% plot(tspan,gamman,lm,'LineWidth',2); 
% 
% figure(3)
% [mutrace,~,mu_o,~,~] = f_carreau(vmaterial,gamma);
% plot(tspan,mutrace/mu_o,lm,'LineWidth',2); 
% 
% figure(4)
% tautrace = mutrace.*gamma;
% stress = po;
% plot(tspan,tautrace/stress,lm,'LineWidth',2); 
% 
% figure(5)
% tautrace = mutrace.*gamma;
% taug = abs([gamma tautrace])';
% taug = sortrows(taug');
% gt = taug(:,1)/fnatural;
% taut = taug(:,2)/stress;
% plot(gt,taut,lm,'LineWidth',2); 