function [X,info]=karcher(A,options)

% X=KARCHER(A1,...,Ap) computes the Karcher mean of positive
%  definite matrices A1,...,Ap, using the relaxed Richardson iteration,
%  where the parameter theta may be chosen automatically,
%  and the initial value is the arithmetic mean
% X=KARCHER(A{1:p}) for a cell-array input
% X=KARCHER(A1,...,Ap,theta) same as above, but theta provided by user
% X=KARCHER(A{1:p},theta) for cell-array input
% Do not use X=KARCHER(A) for a cell-array A, use KARCHER(A{1:p}) instead
%
% varargin: positive definite matrix arguments A1,...,Ap
% theta: the parameter of the iteration
% X: the Karcher mean of A1,...,Ap
% iter: the number of iterations needed by the outer iteration
%
%
% This implementation comes from the Matrix Means toolbox developed by
% D.A. Bini and B. Iannazzo, and available at
% http://bezout.dm.unipi.it/software/mmtoolbox/
%
% References:
% [1] D.A. Bini and B. Iannazzo, "Computing the Karcher mean of symmetric
% positive definite matrices", Linear Algebra Appl., 438-4 (2013),
% pp. 1700-1710.



tic;
p=length(A);
info = struct();

if nargin==1
    options = struct();
end

if ~isfield(options,'MStart')
    options.MStart = sum(cat(3,A{:}),3)/p;
end

if ~isfield(options,'theta')
    options.theta = -1;
end

if ~isfield(options,'stop')
    options.stop = -1;
end

if ~isfield(options,'record')
    options.record = 0;
end

if ~isfield(options,'maxiter')
    options.maxiter = 1000;
end

if options.record
    info.M_rec = cell(1,options.maxiter);
    info.t = zeros(1,options.maxiter);
end


X=options.MStart;
for h=1:p
    R{h}=chol(A{h});
end

if (options.stop == 2)
    tol=1d-15;niold=Inf;
end


for k=1:options.maxiter
    R0=chol(X);
    iR0=inv(R0);
    for h=1:p
        Z=R{h}*iR0;
        [Uz{h} Vz]=schur(Z'*Z);
        V{h}=diag(Vz);
    end
    
    theta = options.theta;
    if (options.theta == -1)
        % automatic choice of theta
        beta=0;gamma=0;
        for h=1:p
            ch=max(V{h})/min(V{h});
            if (abs(ch-1)<0.5)
                dh=log1p(ch-1)/(ch-1);
            else
                dh=log(ch)/(ch-1);
            end
            beta=beta+dh;gamma=gamma+ch*dh;
        end
        theta=2/(gamma+beta);
    end
    
    S=0;
    for h=1:p
        T=Uz{h}*diag(log(V{h}))*Uz{h}';
        S=S+(T+T')/2;
    end
    [Us, Vs]=schur(S);
    Z=diag(exp(diag(Vs*theta/2)))*Us'*R0;
    X=Z'*Z;
    info.M_rec{k} = X;
    info.t(k) = toc;
        
        
%   My stopping criteria
    if (options.stop == 1 && k>=2 && norm(info.M_rec{k}-info.M_rec{k-1},'fro')<= 10^(-12)*norm(info.M_rec{k-1},'fro'))
        %fprintf('Number of iterations Karcher = %d \n',k);
        break;
    end
    

%   Stopping criteria initially proposed in the toolbox    
    if (options.stop == 2)
        ni=max(abs(diag(Vs))); % max(abs(diag(Vs)))=norm(S)
        if ( (ni<norm(X)*tol) || (ni>niold) )
            %fprintf('Number of iterations Karcher = %d \n',k);
            break;
        end
        niold=ni;
    end
    
end

M = info.M_rec{k};
info.tTot = toc;


end