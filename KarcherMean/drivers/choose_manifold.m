function fs = choose_manifold(manifold, retraction, vector_transport)
% This function return some necessary function handles of some manifolds which are required in Riemannian optimization methods.
%
% INPUT:
% retraction : a number represents a retraction. See sphere.m, stiefel.m orthogroup.m and grassmann.m for details
% vector transport : a number represents a vector transport. See sphere.m, stiefel.m orthogroup.m and grassmann.m for details
% OUTPUT:
% Output function handles about the SPD matrix manifold.
%
% Code extracted from http://www.math.fsu.edu/~whuang2/papers/ARLBACMGM.htm
% References: 'A Riemannian Limited-Memory BFGS Algorithm for
% Computing the Matrix Geometric Mean', X. Yuan, W. Huang, P.-A. Absil, K.
% Gallivan.
% (Slight modifications between this version and the original are possible)


fs = SPD(retraction, vector_transport);

end
