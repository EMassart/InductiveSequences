function check_function_handle(fns, name)
% Check whether the required function handle exists or not.
% if it does not exist, an error message is given.
% INPUT:
% fns : a struct that contains function handles.
% name : the name of the required function handle.
% 
% By Wen Huang


% Code extracted from http://www.math.fsu.edu/~whuang2/papers/ARLBACMGM.htm
% References: 'A Riemannian Limited-Memory BFGS Algorithm for
% Computing the Matrix Geometric Mean', X. Yuan, W. Huang, P.-A. Absil, K.
% Gallivan.
% (Slight modifications between this version and the original are possible)

    if(~isfield(fns, name) || ~isa(fcnchk(getfield(fns, name)), 'function_handle'))
        msg = sprintf('Invalid arguments: missing fns.%s or fns.%s is not a function handle', name, name);
        error(msg);
    end
end
