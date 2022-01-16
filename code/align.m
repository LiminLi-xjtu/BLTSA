
function  [mappedX,B,I,obj]= align(X, label,no_dims, k, eig_impl)
%LTSA Runs the local tangent space alignment algorithm
%
%   mappedX = ltsa(X, no_dims, k, eig_impl)
%
% The function runs the local tangent space alignment algorithm on dataset
% X, reducing the data to dimensionality d. The number of neighbors is
% specified by k.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology
for t = 1:length(unique(label))
Branch{t} = find(label==t);
end

    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('eig_impl', 'var')
        eig_impl = 'Matlab';
    end
 
    % Compute neighborhood indices
    if ischar(k)
        warning('Adaptive neighborhood selection often leads to problems in LTSA.');
    end
    disp('Find nearest neighbors...');
    n = size(X, 1);
    [D, ni] = myfind_nn(X, k);
%     [D, maxk,no_dim] = find_nn_adaptive(X);
    for i=1:n
%         I0{i} = ni(i,:);
        Ii = ni(i,:);%intersect(ni(i,:),Branch{label(i)});
%         Ii = find(D(i,:)~=0);
        Ii = Ii(Ii ~= 0);
        if label(i)==2
            Ii = setdiff(Ii,Branch{3});

        end
        if label(i)==3
            Ii = setdiff(Ii,Branch{2});
        end
%         if label(i)==2
%             Ii = setdiff(Ii,union(Branch{3},Branch{4}));
%         end
%         if label(i)==3
%             Ii = setdiff(Ii,union(Branch{2},Branch{4}));
%         end
%         if label(i)==4
%             Ii = setdiff(Ii,union(Branch{2},Branch{3}));
%         end
        I{i} = Ii;
    end


    
    % Compute local information matrix for all datapoints
    disp('Compute local information matrices for all datapoints...');
    Bi = cell(1, n); 
    for i=1:n
        % Compute correlation matrix W
        Ii = I{i};
        kt = numel(Ii);
        Xi = X(Ii,:) - repmat(mean(X(Ii,:), 1), [kt 1]);
        W = Xi * Xi'; 
        W = (W + W') / 2;
        
        % Compute local information by computing d largest eigenvectors of W
        [Vi, Si] = schur(W);
        [s, Ji] = sort(-diag(Si));
		if length(Ji) < no_dims
			no_dims = length(Ji);
			warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
		end
        Vi = Vi(:,Ji(1:no_dims)); 
        
        if 1
%         % Store eigenvectors in G (Vi is the space with the maximum variance, i.e. a good approximation of the tangent space at point Xi)
% 		% The constant 1/sqrt(kt) serves as a centering matrix
		Gi = double([repmat(1 / sqrt(kt), [kt 1]) Vi]);
%         % Compute Bi = I - Gi * Gi'
		Bi{i} = eye(kt) - Gi * Gi';  
        else
            Pi = eye(kt) - Vi*Vi';
            Wi = Pi - repmat(mean(Pi, 1), [kt 1]);
            Bi{i} = Wi*Wi';
        end
    end
    
    % Construct sparse matrix B (= alignment matrix)
    disp('Construct alignment matrix...');
    B = speye(n);
    for i=1:n
%         Ii = intersect(ni(i,:),Branch{label(i)});;
%         Ii = Ii(Ii ~= 0);
        Ii = I{i};
        B(Ii, Ii) = B(Ii, Ii) + Bi{i};							% sum Bi over all points
		B(i, i) = B(i, i) - 1;
    end
	B = (B + B') / 2;											% make sure B is symmetric
	B = full(B);
    
	% For sparse datasets, we might end up with NaNs in M. We just set them to zero for now...
	B(isnan(B)) = 0;
	B(isinf(B)) = 0;
    
    % Perform eigenanalysis of matrix B
	disp('Perform eigenanalysis...');
    tol = 0;
	if strcmp(eig_impl, 'JDQR')
        options.Disp = 0;
        options.LSolver = 'bicgstab';
        [mappedX, D] = jdqr(B, no_dims + 1, tol, options);      % only need bottom (no_dims + 1) eigenvectors
    else
        options.disp = 0;
        options.isreal = 1;
        options.issym = 1;
        [mappedX, D] = eigs(B, no_dims + 1, tol, options);      % only need bottom (no_dims + 1) eigenvectors
    end

    % Sort eigenvalues and eigenvectors
    [D, ind] = sort(diag(D), 'ascend');
    mappedX = mappedX(:,ind);

    % Final embedding coordinates
	if size(mappedX, 2) < no_dims + 1, no_dims = size(mappedX, 2) - 1; end
    mappedX = mappedX(:,2:no_dims + 1);
    
    for i = 1:n
        obj(i) = norm(Bi{i}*mappedX(I{i},:),'fro')^2;
    end
    obj = obj';
%     sumobj = sum(obj);
    