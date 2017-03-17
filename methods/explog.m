function [X,info] = explog(A,options)

% X=EXPLOG(A1,A2,...,AP) computes the explog mean 
%  exp(1/P*(log(A1)+...+log(AP)))
%
% X: the explog mean of A{1},...,A{P}

% Original implementation: The Matrix Mean toolbox, D.A. Bini, B. Iannazzo
% available at http://bezout.dm.unipi.it/software/mmtoolbox/

% Modified by E. Massart

tic;
n=length(A{1});
p=length(A);

S=zeros(n);
for k=1:p
    R = chol(A{k});
    [U, V] = schur(R'*R);
    S=S+U*diag(log(diag(V)))*U';
end

[U, V] = schur(1/p*S);
X = U*diag(exp(diag(V)))*U';
info.tTot = toc; 

end