function [X, F, G, T, timecost, iter,status] = Karcher_mean(As,method,options)

% X = Karcher_mean(As,method,TOL,initial_step) computes the Karcher mean of SPD matrices As1,...,AsK
% INPUT:
%   As:data points
%   method: Riemannian optimization method
%           1 means Riemannian steepest descent method
%           2 means Riemannian conjugate gradient method
%           3 means Riemannian limited-memory BFGS method
%           4 means Riemannian BFGS method
%   TOL: stopping tolerance
%   initial_step: initial iterate
% OUTPUT:
%   X: all iterates
%   F: function values of all iterates
%   G: norm of gradient of all iterates
%   T: accumulated computational time for iterations
%   timecost: total computational time
%   iter: number of iterations

% Code extracted from http://www.math.fsu.edu/~whuang2/papers/ARLBACMGM.htm
% References: 'A Riemannian Limited-Memory BFGS Algorithm for
% Computing the Matrix Geometric Mean', X. Yuan, W. Huang, P.-A. Absil, K.
% Gallivan.
% (Slight modifications between this version and the original are possible)


k = length(As);
n = size(As{1},2);

for i = 1:k
    A.U{i} = As{i};
    A.h{i} = (As{i})^(0.5);
    A.invh{i} = pinv(A.h{i});
    A.inv{i} = pinv(As{i});
end

x1.U = options.MStart;
x1.L = chol(x1.U,'lower');
x1.invL = inv(x1.L);
x1.invU = x1.invL' * x1.invL;
if(~isfield(x1, 'c'))
    for i = 1 : length(As)
        S = A.invh{i} * x1.L;
        [V, D] = schur(S * S');
        x1.c(i) = D(end,end)/D(1,1);
    end
end


d = 0.5 * (n + 1) * n;
H1 = eye(d);
B1 = eye(d);

params.H0 = H1;
params.B0 = B1;
params.x0 = x1;
params.StopCriterion = options.stop;
params.error = options.err;


if(method == 1 || method == 2)
    params.max_t = 120;
else
    params.max_t = options.maxiter;
end


C = max(x1.c);
L = 1 + 0.5 * log(C);
params.initstepsize = 2/(1+L);

params.accuracy = 1e-3;
params.debug = 0;
params.manifold = 0;
params.m = 2;
params.alpha = 1e-4;
params.beta = 0.999;
params.Delta_bar = inf;
params.Delta0 = 1;
params.min_Delta = eps;
params.minstepsize = eps;
params.maxstepsize = 200;
params.ratio = 0.25;
params.restart = 0;
params.retraction = 2;
params.vector_transport = 1;
params.hessian = 1;


if(params.retraction == 1)
    params.linesearch = 2; % Wolfe condition line search
elseif(params.retraction == 2)
    params.linesearch = 1; % Armijo condition line search
end


fns = {};

[X, F, G, T, timecost, iter,status] = driver_SPD(k, A, method, fns, params);


end





