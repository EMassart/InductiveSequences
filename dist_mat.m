function [d,de]=dist_mat(A,B,alg)

% [d,de]=DIST(A,B,alg) computes the distance between A and B in the Riemannian 
%  metric and in the Euclidean metric
%
% A,B: positive definite matrices
% alg: the algorithm for computing the distance
% alg='geigs' uses the generalized eigenvalues of the couple (A,B) 
% alg='chole' uses the algorithm based on the Cholesky factorization
% alg='alter' similar to chole, but using schur decomposition
% d: the Riemannian distance between A and B 
% de: the Euclidean distance between A and B
%
% Note: the algorithm based on the Cholesky factorization is 
% theoretically faster than the one based on the generalized 
% eigenvalues, but in the execution the latter requires less time
% since uses the built-in QZ algorithm

% This implementation comes from the Matrix Means toolbox developed by
% D.A. Bini and B. Iannazzo, and available at
% http://bezout.dm.unipi.it/software/mmtoolbox/

if (nargin<3)
  alg='geigs';
end

if (strcmp(alg,'alter'))
  mA=1/rcond(A);mB=1/rcond(B);
  if (mA>mB) % swap A and B if B is better conditioned
    C=A;A=B;B=C;t=1-t;
  end	
  RA=chol(A);RB=chol(B);

  Z=RB/RA;
  [U V]=schur(Z'*Z);
  d=norm(log(diag(V)));
end

if (strcmp(alg,'chole'))
  mA=1/rcond(A);mB=1/rcond(B);
  if (mA>mB) % swap A and B if B is better conditioned
    C=A;A=B;B=C;
  end	
  RA=chol(A);RB=chol(B);
  Z=RB/RA;
  d=norm(log(eig(Z'*Z)));
end

if (strcmp(alg,'geigs'))
  d=norm(log(eig(A,B)));
end

if (nargout>1)
  de=norm(A-B,'fro');
end