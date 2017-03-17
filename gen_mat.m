function [ A ] = gen_mat( problem )
% gen_mat(problem) : generates a set of SPD matrices 
% of given size, with norm equal to one.
%
% problem is a structure containing the fields : 
%    problem.number: number of data matrices
%    problem.size: size of the matrices
%    (problem.cond): magnitude order for the logarithm of the condition number of the
%                    data. If problem.cond is empty, then the matrices are
%                    generated according to the Wishart definition.
%    (problem.version): describes how to generate SPD matrices with 
%                       approximate) given condition number.

% Author: E. Massart


A = cell(1,problem.number);

for i_mat = 1:problem.number
    
if (~isfield(problem,'cond'))           % matrices generated according to Wishart distribution
    A{i_mat} = randn(problem.size,problem.size);
    A{i_mat} = A{i_mat}'*A{i_mat};
else
    
    if ~isfield(problem, 'version')
        problem.version = 3;
    end
    
    [Q,~] = qr(randn(problem.size));

    % First variant: one eigenvalue considerably smaller than the others
    if problem.version == 1
        D = diag([rand(1,problem.size-1)+1,10^(-problem.cond)]);
        A{i_mat} = Q*D*Q';
        
    % Second variant: one eigenvalue considerably bigger than the others
    elseif problem.version == 2
        D = diag([1,10^(-problem.cond)+10^(-problem.cond)*rand(1,problem.size-1)]);
        A{i_mat} = Q*D*Q';
        
    % Third variant: two sets of approximatively the same number of eigenvalues, the eigenvalues 
    % of the first set are problem.cond orders of magnitude higher than those of the second set.
    elseif problem.version == 3
        D = diag([1+rand(1,floor(problem.size/2)),10^(-problem.cond)*(1+rand(1,problem.size - floor(problem.size/2)))]);
        A{i_mat} = Q*D*Q';
    end
end

nor = norm(A{i_mat},'fro');
A{i_mat} = A{i_mat}./nor;

lambdas = eig(A{i_mat});
e = find(lambdas<=0);

if(size(e)~= 0)
    disp('The matrix is not positive definite');
end

end

end



