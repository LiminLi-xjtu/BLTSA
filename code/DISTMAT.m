%X: m*n
%Dis:  m*m

function Dis = DISMAT(X)
[m,n] = size(X);
A = X*X';
l = diag(A);
Dis = sqrt(repmat(l,[1,m])+ repmat(l',[m,1]) - 2*A);
							  
end
