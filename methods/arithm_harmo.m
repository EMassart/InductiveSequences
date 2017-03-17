function [M,info] = arithm_harmo(As,options)
%Computes the Arithmetic-harmonic mean of the matrices contained in the
%cell array As

% Author: E. Massart

tic;

k = length(As);

Ar = As{1};
Ha = inv(As{1});
for i = 2:k
    Ar = Ar + As{i};
    Ha = Ha + inv(As{i});
end
Ar = Ar/k;
Ha = inv(Ha/k);

% c_M = [cond(Ar), cond(Ha)];
    
M = sharp(Ar,Ha,0.5);
info.tTot = toc;  
end