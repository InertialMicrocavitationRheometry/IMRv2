figure(1)
tspan = Tout*omegan;
Roft = (Rout(:,1)/Ro_w)-1;
plot(tspan,Roft, lm,'LineWidth',2); 

figure(2)
gamma = -4*(Rout(:,2))./(Rout(:,1));
gamman = gamma;
plot(tspan,gamman,lm,'LineWidth',2); 

figure(3)
mutrace = carreau(vmaterial,gamma).*ones(size(gamma));
plot(tspan,mutrace/mu_o,lm,'LineWidth',2); 

figure(4)
tautrace = mutrace.*gamma;
plot(tspan,tautrace,lm,'LineWidth',2); 

figure(5)
tautrace = mutrace.*gamma;
taug = abs([gamma tautrace])';
taug = sortrows(taug');
gt = taug(:,1);
taut = taug(:,2);
plot(gt,taut,lm,'LineWidth',2); 