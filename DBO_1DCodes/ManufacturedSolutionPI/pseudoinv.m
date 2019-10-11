function [NewMat,MatInv, NewY] = pseudoinv(Mat, Y)

global tol_PI
[U, S, V] = svd(Mat);
N = size(Mat,1);
Sinv = zeros(N,N);
% Convert Y_DO to Y_BO 
YBO = inv(sqrt(S))*V'*Y'; 
for i=1:N
   if(S(i,i) < tol_PI)
       S(i,i) = tol_PI;
       Sinv(i,i) = 1/tol_PI;
   else
       S(i,i)=S(i,i);
       Sinv(i,i) = 1/S(i,i);
   end
end
NewMat = U*S*V';
MatInv = U*Sinv*V';

NewY = V*sqrt(S)*YBO;
NewY = NewY';
