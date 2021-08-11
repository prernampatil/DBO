% Plot the error in singular values for DO and DBO
clc
clear 

load('Dirichlet_d8_r9.mat')

ErrDBO = sqrt((SDNS-Sigma).^2);
ErrDO  = sqrt((SDNS-SDO).^2);

figure(12)

plot(T(1:500:end),ErrDBO(3,1:500:end),'-.v','color',RGB1,'MarkerFaceColor',RGB1,LW,1.5);
hold on
plot(T(1:500:end),ErrDO(3,1:500:end),'-.d','color',RGB2,'MarkerFaceColor',RGB2,LW,1.5);

plot(T(1:500:end),ErrDBO(6,1:500:end),'-.o','color',RGB1,'MarkerFaceColor',RGB1,LW,1.5);
hold on
plot(T(1:500:end),ErrDO(6,1:500:end),'-.*','color',RGB2,'MarkerFaceColor',RGB2,LW,1.5);

plot(T(1:500:end),ErrDBO(9,1:500:end),'--','color',RGB1,LW,1.5);
hold on
plot(T(1:500:end),ErrDO(9,1:500:end),':','color',RGB2,LW,2.5);

legend('DBO \Sigma_{i=3}','DO    \Sigma_{i=3}','DBO \Sigma_{i=6}','DO    \Sigma_{i=6}',...
       'DBO \Sigma_{i=9}','DO    \Sigma_{i=9}','Location','northeast','NumColumns',3)
set(gca,'Yscale','log');
ylim([10^-14 10^1])
xlabel('Time');
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
ylabel('Singular Values Error');
set(gca,'FontName','Times New Roman','FontSize',18);
grid on
ax= gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'off';