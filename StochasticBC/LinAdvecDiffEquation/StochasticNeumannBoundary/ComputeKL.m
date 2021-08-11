function [U,S,V] = ComputeKL(A,wp,wr)
global NModes 
C = A'*diag(wp)*A; C= (C+C')/2;
[V,L] = eig(C*diag(wr)); V = real(V); L = real(L);
L     = diag(L);
[L, I] = sort(L,'descend');
V = V(:,I(1:NModes));
C = abs(diag(diag(V'*diag(wr)*V)));
V = V/sqrt(C);
S= diag(real(sqrt(L(1:NModes))));
U = A*diag(wr)*V/S;