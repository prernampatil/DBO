% Plot the singular values 
% Written by: Prerna Patil 

clc
clear 
close all 

% Singular values 
load('BurgersDBC_d4_r4.mat');
Tf = 5;

figure(1)
for i = 1:NModes
    semilogy(dt*(1:nTimeStep/nSkip)*nSkip,SDNS(i,:),'--k',LW,1.5);
    hold on 
    semilogy(dt*(1:nTimeStep),SDO(i,:),'--','color',RGB2,LW,1.5);
    semilogy(dt*(1:500:nTimeStep),SDO(i,1:500:nTimeStep),'*','color',RGB2,LW,1.5);
    
    semilogy(dt*(1:nTimeStep),Sigma(i,:),':','color',RGB1,LW,1.5);
    semilogy(dt*(250:500:nTimeStep),Sigma(i,250:500:nTimeStep),'+','color',RGB1,LW,1.5);
    
end
xlabel('Time')
ylabel('Singular Values');


% Singular values 
clear
load('BurgersDBC_d4_r6.mat');
Tf = 5;

figure(2)
for i = 1:NModes
    semilogy(dt*(1:nTimeStep/nSkip)*nSkip,SDNS(i,:),'--k',LW,1.5);
    hold on 
    semilogy(dt*(1:nTimeStep),SDO(i,:),'--','color',RGB2,LW,1.5);
    semilogy(dt*(1:500:nTimeStep),SDO(i,1:500:nTimeStep),'*','color',RGB2,LW,1.5);
    
    semilogy(dt*(1:nTimeStep),Sigma(i,:),':','color',RGB1,LW,1.5);
    semilogy(dt*(250:500:nTimeStep),Sigma(i,250:500:nTimeStep),'+','color',RGB1,LW,1.5);
    
end
xlabel('Time')
ylabel('Singular Values');


% Singular values 
clear
load('BurgersDBC_d4_r8.mat');
Tf = 5;

figure(3)
for i = 1:NModes
    semilogy(dt*(1:nTimeStep/nSkip)*nSkip,SDNS(i,:),'--k',LW,1.5);
    hold on 
    semilogy(dt*(1:nTimeStep),SDO(i,:),'--','color',RGB2,LW,1.5);
    semilogy(dt*(1:500:nTimeStep),SDO(i,1:500:nTimeStep),'*','color',RGB2,LW,1.5);
    
    semilogy(dt*(1:nTimeStep),Sigma(i,:),':','color',RGB1,LW,1.5);
    semilogy(dt*(250:500:nTimeStep),Sigma(i,250:500:nTimeStep),'+','color',RGB1,LW,1.5);
    
end
xlabel('Time')
ylabel('Singular Values');
h= zeros(3,1);
h(1) = plot(NaN,NaN,'--k','linewidth',1.5);
h(2) = plot(NaN,NaN,':+','color',RGB1,LW,2.5);
h(3) = plot(NaN,NaN,'-*','color',RGB2,LW,1.5);
legend(h,'KL','DBO','DO','Location','southeast' );
set(gca,'FontName','Times New Roman','FontSize',18);
box on 
grid on 
ax= gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';

figure(1)
hold on 
for i = 6:NModes
    semilogy(dt*(1:nTimeStep/nSkip)*nSkip,SDNS(i,:),'--k',LW,1.5);
    hold on 
end
ylim([10^-8 10^1]);
xlim([Ts Tf]);
h= zeros(3,1);
h(1) = plot(NaN,NaN,'--k','linewidth',1.5);
h(2) = plot(NaN,NaN,':+','color',RGB1,LW,2.5);
h(3) = plot(NaN,NaN,'-*','color',RGB2,LW,1.5);
legend(h,'KL','DBO','DO' ,'Location','southeast' );
set(gca,'FontName','Times New Roman','FontSize',18);
box on 
grid on 
ax= gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));


figure(2)
hold on 
for i = 8:NModes
    semilogy(dt*(1:nTimeStep/nSkip)*nSkip,SDNS(i,:),'--k',LW,1.5);
    hold on 
end
ylim([10^-8 10^1]);
xlim([Ts Tf]);
h= zeros(3,1);
h(1) = plot(NaN,NaN,'--k','linewidth',1.5);
h(2) = plot(NaN,NaN,':+','color',RGB1,LW,2.5);
h(3) = plot(NaN,NaN,'-*','color',RGB2,LW,1.5);
legend(h,'KL','DBO','DO','Location','southeast'  );
set(gca,'FontName','Times New Roman','FontSize',18);
box on 
grid on 
ax= gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));

figure(3)
ylim([10^-8 10^1]);
xlim([Ts Tf]);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));