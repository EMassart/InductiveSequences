function [M,info] = arithm(As,options)
%Computes the Arithmetic mean of the matrices contained in the
%cell array As

% Author: E. Massart

tic;
M = sum(cat(3,As{:}),3)/length(As);
info.tTot = toc;   
end