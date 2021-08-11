NModes=7;
for i=1:NModes
    figallM = figure(2);
    figallM = subplot(9,1,i);
    A = figallM.Position;
    if(i==1)
        legA = legend;
        legall = legA.Position;
        
    end
    figchange = figure(1);
    figchange = subplot(NModes,1,i);
    figchange.Position = A;
    if(i==1)
        leg = legend;
        leg.Position= legall;
    end
end
