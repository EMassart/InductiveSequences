function [ M, info ] = Inductive( A, options)
%%Compute the Inductive sequence of the data points A{1}, ..., A{N}.
% Options (default values are indicated in []):
%   option.maxiter [50]: number of times that the algorithm visits each data
%                        point.
%   option.record [0]: 
%      = 1 to record the iterates and elapsed time per iteration
%      = 0 otherwise
%   option.shuffling ['shuffling']: 
%      shuffling algorithm used to  generate the sequence of permutations
%      = 'shuffling': implementation of Algorithm 1 of the paper "Matrix
%                     geometric means based on shuffled inductive sequences", 
%                     E.Massart, J. Hendrickx, P.-A. Absil
%      = 'no': repeat the sequence 1,2,...,N
%      = 'random': succession of random permutations of 1,2,...,N
%      = 'random_replacement': succession of indices taken with replacement
%        among 1,...,N (the same index might be taken several times
%        consecutively)

% Author: E. Massart

tic;
N = length(A);
info = struct();

if nargin==1
    options = struct();
end

if ~isfield(options,'record')
    options.record = 0;
end

if ~isfield(options,'maxiter')
    options.maxiter = 50;
end

if ~isfield(options,'shuffling')
    options.shuffling = 'shuffling';
end

if options.record  
    info.M_rec = cell(1,options.maxiter);
    info.t = zeros(1,options.maxiter);
end

% Factorization of the data
% To reduce computations, we directly work on the square roots
% (Cholesky factorizations) of the matrices
R = cell(1,N);
R_inv = cell(1,N);
for i = 1:N
    R{i} = chol(A{i});
    R_inv{i} = inv(R{i});
end  

%Initial point
RM = R{1};
options_Shuff.count = 1;
    

for i = 1:options.maxiter
    
    switch options.shuffling
        case 'no'
            G = 1:N;
        case 'random'
            G = randperm(N);
        case 'random_replacement'
            G = randi(N,1,N);
        case 'shuffling'
            if (i==1)
                options_Shuff.count = 1;
                G = 1:N;
            elseif mod(i,2)==1
                G = G(end:-1:1);
                [G,options_Shuff] = riffle_shuffle(G,options_Shuff);
            else
                G = G(end:-1:1);
            end
    end
    
    % Makes N steps along geodesics, according to permutation G.
    for i_loc = 1:N
        indx = N*(i-1)+i_loc;           % steps counter (used to for step length computation)
        Z = RM*R_inv{G(i_loc)};         
        [U V] = schur(Z'*Z);
        RM=diag(diag(V).^((indx-1)/(2*indx)))*U'*R{G(i_loc)};       
    end

    % Records the data
    if options.record
        info.M_rec{i} = RM;
        info.t(i) = toc;
    end
    
end


% Final result
M = RM'*RM;
info.tTot = toc;


% Computes the intermediary matrices from the square roots that have been
% previously recorded.
if options.record
    for i = 1:length(info.M_rec)
        info.M_rec{i} = info.M_rec{i}'*info.M_rec{i};
    end
end


end