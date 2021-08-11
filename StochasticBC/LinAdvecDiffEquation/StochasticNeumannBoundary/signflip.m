for i=1:NModes
    

    % Change sign for K
    if(i>1)
        
        tempvec = find(diff(sign(BCModesDBO(i,:))))+1;
        p =1;
        while(p<length(tempvec))
            if(abs(BCModesDBO(i,tempvec(p))-BCModesDBO(i,tempvec(p)-1))>0.1)
                BCModesDBO(i,tempvec(p):tempvec(p+1)-1) = -BCModesDBO(i,tempvec(p):tempvec(p+1)-1);
            end
            p = p+1;
        end
    end
    sgnKL = sign(BCModesKL(i,:));
    sgnDBO = sign(BCModesDBO(i,:));
    signchange = -0.5*abs(sgnKL-sgnDBO);
    signchange(signchange==0) = 1;
    BCModesKL(i,:) = BCModesKL(i,:).*signchange;
    subplot(NModes,1,i)
    plot(T(1:50:end), BCModesKL(i,1:50:end),'-','color','black',LW,1.5);
    hold on
    plot(T(1:50:end), BCModesDBO(i,1:50:end),':','color', RGB1, LW ,2.5)
    xticks([]);
    
    if(i==NModes)
        xlabel('Time')
        xticks([1 2 3 4 5 6 7 8 9 ])
        xticklabels({'1','2','3','4','5','6','7','8','9'})
        xlim([Ts Tf])
    end
    ylim([-2 2]);
    ylabel(sprintf('$u_{%d}(0,t)$', i))
    if(i==1)
        legend('KL','DBO');
    end
    set(gca,'FontName','Times New Roman','FontSize',15);
end