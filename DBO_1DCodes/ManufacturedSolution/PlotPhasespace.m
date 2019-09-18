function PlotPhasespace(DO, BO, DBO, YDO, YBO, YDBO, YExact)
clf
global RGB1 RGB2 RGB3 
LW = 'linewidth';
% plot the Y phase space for the solutions
if(DO)
    plot(YDO(1:8:end,1), YDO(1:8:end,2), '+','color', RGB3, LW, 1.5);
    hold on
end
if(BO)
    plot(YBO(1:8:end,1), YBO(1:8:end,2), 'o', 'color', RGB1, LW, 1.5);
    hold on
end
if(DBO)
    plot(YDBO(1:8:end,1), YDBO(1:8:end,2), 'x', 'color', RGB2, LW, 1.5);
    hold on
end
plot(YExact(1:8:end,1), YExact(1:8:end,2), 'sqk', LW, 1.5);
xlabel('$\mathrm{Y_1}$', 'interpreter','latex');
ylabel('$\mathrm{Y_2}$', 'interpreter','latex');
xticks([-1.5  -0.5  0.5  1.5])
xticklabels({'-1.5','-0.5', '0.5',  '1.5'})
yticks([-2.5  -1.5 -0.5  0.5  1.5 2.5])
ylim([-2.5 2.5])
xlim([-1.5 1.5])
yticklabels({'-2.5', '-1.5','-0.5', '0.5',  '1.5', '2.5'})
title('Phase space Y_1 vs. Y_2', 'Interpreter', 'tex');
h= zeros(4,1);
h(1) = plot(NaN, NaN, 'o','color', RGB1, LW, 1.5);
h(2) = plot(NaN, NaN, '+','color', RGB3, LW, 1.5);
h(3) = plot(NaN, NaN, 'x','color', RGB2, LW, 1.5);
h(4) = plot(NaN, NaN, 'sqk', LW, 1.5);
legend(h,'BO', 'DO', 'DBO', 'Analytical');
set(gca, 'FontSize', 20, 'Fontname', 'Times New Roman');
end