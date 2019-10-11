function PlotSpatialModes(DO,PIDO, BO, DBO, Lambda, EigVec, UBO, UDO, UPIDO, UDBO, UExact, M )
global RGB1 RGB2 RGB3 RGB4
global epname
global x 
LW = 'linewidth';
if(BO)
    UBOtoDO = (inv(sqrt(Lambda))*EigVec'*UBO')';
    plot(x, UBOtoDO(:,M), ':','color', RGB1, LW, 2.5);
    hold on
end
if(DO)
    plot(x(1:20:end), UDO(1:20:end,M), 'd','color', RGB3,'MarkerSize',3.5);
    hold on
    plot(x, UDO(:,M), '-','color', RGB3, LW, 1.5);
end
if(PIDO)
    plot(x(10:20:end), UPIDO(10:20:end,M), 'p','color', RGB4,'MarkerSize',5.5);
    hold on
    plot(x, UPIDO(:,M), '-','color', RGB4, LW, 1.5);
end
if(DBO)
    plot(x, UDBO(:,M),'-','color', RGB2, LW, 1.5);
    hold on
    plot(x(5:20:end), UDBO(5:20:end,M),'x','color', RGB2, LW, 1.5);
end
plot(x, UExact(:,M), '.k', LW, 1.5)
xlim([0, 2*pi])
ylim([-1.2 1.2])
xticks([0  pi/2  pi  3*pi/2  2*pi])
xticklabels({'0','\pi/2', '\pi',  '3\pi/2',  '2\pi'})
xlabel('X')
h= zeros(5,1);
h(1) = plot(NaN, NaN, '--k',LW, 1.5);
h(2) = plot(NaN, NaN, ':','color', RGB1, LW, 2.5);
h(3) = plot(NaN, NaN, '-d','color', RGB3, LW, 1.5,'MarkerSize',3.5);
h(4) = plot(NaN, NaN, '-p','color', RGB4, LW, 1.5,'MarkerSize',5.5);
h(5) = plot(NaN, NaN, '-x','color', RGB2, LW, 1.5);
legend(h,'Analytical','BO', 'DO', 'PI-DO', 'DBO');
set(gca, 'FontSize', 20, 'Fontname', 'Times New Roman')

end