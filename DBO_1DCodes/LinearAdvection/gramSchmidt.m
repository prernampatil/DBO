function [Y, U] = gramSchmidt(bY, bU)
global wr wp
[~,N] = size(bU);
Y(:,1) = bY(:,1)/sum(bY(:,1).*bY(:,1).*wr);
U(:,1) = bU(:,1)/sum(bU(:,1).*bU(:,1).*wp);
for i=2:N
    tempY = bY(:,i);
    tempU = bU(:,i);
    for j=1:i-1
        tempY = tempY - sum( bY(:,i).*bY(:,j).*wr )/sum( bY(:,j).*bY(:,j).*wr )*bY(:,j);
        tempU = tempU - sum( bU(:,i).*bU(:,j).*wp )/sum( bU(:,j).*bU(:,j).*wp )*bU(:,j);
    end
    tempY = tempY/ sum(tempY.*tempY.*wr);
    Y(:,i) = tempY;
    tempU = tempU/ sum(tempU.*tempU.*wp);
    U(:,i) = tempU;
end