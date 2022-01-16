addpath('data')
addpath('code')

%use Dimensionality.R to reduce the dimensionality

load sim
root = 1;
k = 200;
epi = 0.35;
label1 = mylabelling(X, 1, k,root,epi);
[Y1,B1,I1,obj1] = myalign(X,label1, 1, k);
Y1 = Y1/(max(Y1)-min(Y1));
Y1 = Y1+abs(min(Y1));
scatter(X(:,1),X(:,2),25,Y1),hold on
title('BLTSA')%标题