
function  [label,S] = mylabelling(X, no_dims, k,root,epi)

% 
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
    [n,d] = size(X);
    [D, ni] = myfind_nn(X, k);
    for i=1:n       
        Ii = ni(i,:);    
        Ii = Ii(Ii ~= 0);
        I{i} = Ii;
    end
    
    for i=1:n
        % Compute correlation matrix W
        Ii = I{i};
        kt = numel(Ii);
        Xi = X(Ii,:) - repmat(mean(X(Ii,:), 1), [kt 1]);
        Xii = X(Ii,:) - repmat(X(i,:), [kt 1]);
        W = Xii * Xii'; 
%         W = (W + W') / 2;
        a(i) = length(find(W(:)<0));  
        % Compute local information by computing d largest eigenvectors of W
        [Ui,Si,Vi] = svds(Xi',no_dims+1,'largest');
        U{i} = Ui(:,1:no_dims);
        si = diag(Si);
        Z{i} = Si(1:no_dims,1:no_dims)*Vi(:,1:no_dims)';        
        sigma(i) = si(end)/sum(si);
    end


    % determine end points and label them
    End = find(a==0); 
    End = [root setdiff(End,root)];   

    num_branch = length(End);

    % determine nonbranch points and label them
    NonBranch = find(sigma<epi);    
    bp =  find(sigma>=epi);    
%     figure,plot(X(:,1),X(:,2),'g*'),hold on
%     plot(X(bp,1),X(bp,2),'r*')
    DD = myfind_nn(X(NonBranch,:), 40);    
    DD = sparse(DD);
    cidx = sc(DD,median(DD(:)),num_branch);
%     cidx = sc(DD,0,num_branch);
    b(NonBranch) = cidx;
    
    
    
    
    for t = 1:num_branch
        S{t} = find(b==b(End(t)));
        label(S{t}) = t;
    end
    Branch = setdiff(1:n,NonBranch);
    

    % extend the labelling to branch points
    D0 = DISTMAT(X);
    bb = [];
    while length(Branch)>0
    

%     step 1: find the easiest branch point


    t = randperm(3); t = t(1);
    Dt = D0(Branch,S{t});
    [a,b] = min(Dt(:));
    [bid,cid] = ind2sub(size(Dt),b);
    branch = Branch(bid);
    
    
    
    %step2: classify current_branch
    % by tangent space
    for t = 1:num_branch
       dt = D0(branch,S{t});
       [~,b] = min(dt);
       near(t) = S{t}(b);
       T = U{near(t)};
       X(branch,:);
       dist(t) = norm((eye(d)-T*T')*(X(branch,:)'-X(near(t),:)'),2);
    end    
%     
%     plot(X(near(1),1),X(near(1),2),'bo'),hold on
%     plot(X(near(2),1),X(near(2),2),'ro'),hold on
%     plot(X(near(3),1),X(near(3),2),'go'),hold on
%     for t = 1:num_branch
%         near(t)
%         text(X(near(t),1),X(near(t),2),num2str(near(t)),'FontSize',20);
%     end
%     hold on;
    
    [~,bb] = min(dist);
    label(branch) = bb;
    S{bb} = [S{bb} branch];
%     if b == 1
%         plot(X(branch,1),X(branch,2),'bo'); 
%     elseif b==2
%         plot(X(branch,1),X(branch,2),'ro');
%     elseif b==3
%         plot(X(branch,1),X(branch,2),'go');
%     end
%     hold off
    
    % update branch and nonbranch 
    Branch = setdiff(Branch,branch);
    NonBranch = [NonBranch branch];
    %update tagent space
    I{branch} = intersect(I{branch},S{bb});
    Ii = I{branch}; kt = numel(Ii);
    Xi = X(Ii,:) - repmat(mean(X(Ii,:), 1), [kt 1]);
   % Compute local information by computing d largest eigenvectors of W
   [Ui,Si,Vi] = svds(Xi',no_dims+1,'largest');
   U{i} = Ui(:,1:no_dims);
    
    end 

end



