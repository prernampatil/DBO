function EYYY = Compute_EYYY(Y)

global wr
global NModes

for i=1:NModes
    for j=1:NModes
        for k=1:NModes 
            EYYY(i,j,k) = wr'*(Y(:,i).*Y(:,j).*Y(:,k));
        end
    end
end
