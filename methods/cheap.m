function [X,info]=cheap(A,options)

% [X,iter]=CHEAP(A1,A2,...,AP) computes the mean defined by Bini and Iannazzo
%  in [1]
% [X,iter]=CHEAP(A{1:p}) for a cell-array input
%
% X: the Cheap mean of A{1},...,A{p}
% 
% References
% [1] D.A. Bini and B. Iannazzo, "A note on computing matrix geometric 
% means", Adv. Comput. Math., 35-2/4 (2011), pp. 175-192.

% Original implementation: The Matrix Mean toolbox, D.A. Bini, B. Iannazzo
% available at http://bezout.dm.unipi.it/software/mmtoolbox/

% Modified by E. Massart

tic;
n=length(A{1});
p=length(A);

if nargin==1
    options = struct();
end

if ~isfield(options,'record')
    options.record = 0;
end

if ~isfield(options,'maxiter')
    options.maxiter = 10;
end

if options.record
    info.M_rec = cell(1,options.maxiter);
    info.t = zeros(1,options.maxiter);
end


for k=1:options.maxiter
  for h=1:p
    R{h}=chol(A{h});
    RI{h}=inv(R{h});
  end
  for h=1:p
    S=zeros(n);
    for ell=1:p
      if (ell~=h)
        % computes S=S+logm(RI{h}'*A{ell}*RI{h})
        Z=R{ell}*RI{h};
        [U V]=schur(Z'*Z);
        T=U*diag(log(diag(V)))*U';
        S=S+(T+T')/2;
      end
    end
    [U V]=schur(1/p*S);
    T=diag(exp(diag(V))).^(1/2)*U'*R{h};
    A1{h}=T'*T;
  end
    
  for h=1:p
    A{h}=A1{h};
  end
    
  info.M_rec{k} = A1{1};
  for h=2:p
      info.M_rec{k} = info.M_rec{k}+A1{h};
  end
  info.M_rec{k} = info.M_rec{k}/p;
  info.t(k) = toc;
end

X = info.M_rec{k};
info.tTot = toc;

end