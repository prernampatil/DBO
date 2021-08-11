% Code to post process the error plots of the Boundary error 
% Written by: Prerna Patil 
% Date: 12th May 2021 

clc
clear 

% Boundary Errors 
load('BurgersDBC_d4_r4.mat');
Tf =5; 

% Plot for r=5 
figure(5)
hold on
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndry(nTswitch+1:nTimeStep-3),'-','color',RGB1,LW,1.5);
semilogy(dt*(nTswitch+500:1000:nTimeStep-3),ErrBndry(nTswitch+500:1000:nTimeStep-3),'x','color',RGB1,...
                                            'MarkerFaceColor',RGB1,'MarkerSize',8,LW,2.0);

hold on 
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndryDO(nTswitch+1:nTimeStep-3),':','color',RGB2,LW,2.5);
semilogy(dt*(nTswitch+1:1000:nTimeStep-3),ErrBndryDO(nTswitch+1:1000:nTimeStep-3),'v','color',RGB2,'MarkerFaceColor',RGB2,LW,1.5);

clear
load('BurgersDBC_d4_r6.mat');
Tf =5; 
figure(5)
hold on
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndry(nTswitch+1:nTimeStep-3),'-','color',RGB1,LW,1.5);
semilogy(dt*(nTswitch+500:1000:nTimeStep-3),ErrBndry(nTswitch+500:1000:nTimeStep-3),'o','color',RGB1,'MarkerFaceColor',RGB1,LW,1.0);

hold on 
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndryDO(nTswitch+1:nTimeStep-3),':','color',RGB2,LW,2.5);
semilogy(dt*(nTswitch+1:1000:nTimeStep-3),ErrBndryDO(nTswitch+1:1000:nTimeStep-3),'*','color',RGB2,'MarkerFaceColor',RGB2,'MarkerSize',6,...
                            LW,1.5);

clear
load('BurgersDBC_d4_r8.mat');
Tf =5; 
figure(5)
hold on
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndry(nTswitch+1:nTimeStep-3),'-','color',RGB1,LW,1.5);
semilogy(dt*(nTswitch+500:1000:nTimeStep-3),ErrBndry(nTswitch+500:1000:nTimeStep-3),'d','color',RGB1,'MarkerFaceColor',RGB1,LW,1.5);

hold on 
semilogy(dt*(nTswitch+1:nTimeStep-3),ErrBndryDO(nTswitch+1:nTimeStep-3),':','color',RGB2,LW,2.5);
semilogy(dt*(nTswitch+1:1000:nTimeStep-3),ErrBndryDO(nTswitch+1:1000:nTimeStep-3),'x','color',RGB2,'MarkerSize',8,...
                    'MarkerFaceColor',RGB2,LW,2.5);

xlabel('Time');
ylabel('$\mathcal{E}_b$');
set(gca,'FontName','Times New Roman','FontSize',18);
ylim([10^-6 10^0])
xlim([Ts Tf]);
box on 
set(gca,'Yscale','log');
ax= gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
grid on 
h = zeros(6,1);
h(1) = semilogy(NaN,NaN,'-x','color',RGB1,'MarkerFaceColor',RGB1,'MarkerSize',8,LW,1.5);
h(2) = semilogy(NaN,NaN,':v','color',RGB2,'MarkerFaceColor',RGB2,'MarkerSize',6,LW,2.5);
h(3) = semilogy(NaN,NaN,'-o','color',RGB1,'MarkerFaceColor',RGB1,'MarkerSize',6, LW,1.5);
h(4) = semilogy(NaN,NaN,':*','color',RGB2,'MarkerFaceColor',RGB2,'MarkerSize',6,LW,2.5);
h(5) = semilogy(NaN,NaN,'-d','color',RGB1,'MarkerFaceColor',RGB1,'MarkerSize',6,LW,1.5);
h(6) = semilogy(NaN,NaN,':x','color',RGB2,'MarkerFaceColor',RGB2,'MarkerSize',8,LW,2.5);
legend(h,{'DBO r=4','DO    r=4','DBO r=6','DO    r=6','DBO r=8','DO    r=8'},'Location','southeast');
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
