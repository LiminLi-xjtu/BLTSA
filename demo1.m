addpath('data')
addpath('code')

%use Dimensionality.R to reduce the dimensionality

X = csvread('sim.csv',1,1);                %data after dimensionality reduction
root = 1;
k = 35;
epi = 0.45;
label = mylabelling(X, 1, k,root,epi);      %calculate labels
[Y1,B1,I1,obj1] = align(X,label, 1, 180);   %calculate pseudotime
Y1 = Y1/(max(Y1)-min(Y1));
Y1 = Y1+abs(min(Y1));
scatter(-X(:,1),X(:,2),25,Y1),hold on