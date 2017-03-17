
% Basic example: generates a set of data points and compute their mean using
% the implementations proposed 

% Author: E. Massart

clear all;
close all; 
clc;

problem.size = 10;
problem.number = 10;
problem.cond = 8;

A = gen_mat(problem);
options = struct();

M{1} = arithm(A,options);
M{2} = arithm_harmo(A,options);
M{3} = cheap(A,options);
M{4} = explog(A,options);
M{5} = Inductive(A,options);
M{6} = Inductive_parallel(A,options);
M{7} = LRBFGS_toolbox(A,options);

M{8} = karcher(A,options);

for i = 1:7
    d(i) = dist_mat(M{i},M{8});
end