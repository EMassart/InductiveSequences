function [X_out,info] = LRBFGS_toolbox( A, options )
% Runs the LRBFGS algorithm to compute the Karcher mean of data points
% A{1},...,A{N}.
%
% This function calls the 'Karcher_mean' function which was initially 
% proposed in http://www.math.fsu.edu/~whuang2/papers/ARLBACMGM.htm
% References:  'A Riemannian Limited-Memory BFGS Algorithm for
% Computing the Matrix Geometric Mean', X. Yuan, W. Huang, P.-A. Absil, K.
% Gallivan.

% Updated on 26/10/2016:
%     - added one output to the function
%     - added default values for fields of input "options".
%
% Author: E. Massart

if nargin==1
    options = struct();
end

if ~isfield(options,'maxiter')
    options.maxiter = 50;
end

if ~isfield(options,'MStart')
    options.MStart = sum(cat(3,A{:}),3)/length(A);
end

if ~isfield(options,'err')
    options.err = 10^(-12);
end

if ~isfield(options,'stop')
    options.stop = 5;
end


[X, F, G, T, timecost, iter,status] = Karcher_mean(A,3,options);

info.M_rec = X;
info.t = T;                     
info.status = status;          
X_out = X{end}.U;

end