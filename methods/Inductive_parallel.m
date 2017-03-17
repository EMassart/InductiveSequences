function [ M, info ] = Inductive_parallel( A, options)
% Parallel variant of Inductive sequences.
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
%   options.mean ['Arith']:
%      mean used to perform the ast averaging step
%      = 'Arith': arithmetic mean
%      = 'AH': arithmetic-harmonic mean
%      = 'Inductive': Inductive mean

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

if ~isfield(options,'mean')
    options.mean = 'Arith';
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

if options.record  
    info.M_rec = cell(1,options.maxiter);
    info.t = zeros(1,options.maxiter);
end

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
    
    % Computes a new matrix B_i, i.e., the Inductive mean of A{1},...,A{N}
    % ordered according to G.
    RM = R{G(1)};
    for i_loc = 2:N
        indx = i_loc;
        Z = RM*R_inv{G(i_loc)};
        [U, V] = schur(Z'*Z);
        RM=diag(diag(V).^((indx-1)/(2*indx)))*U'*R{G(i_loc)};
    end
    B = RM'*RM;
    
    % Updates the mean of the matrices B_i
    if i==1
        M = B;
        if strcmp(options.mean,'AH')
            Ma = B;
            Ma_inv = inv(B);
        end
    elseif strcmp(options.mean,'Arith')
        M = (i-1)*M/i+ B/i;
    elseif strcmp(options.mean,'AH')
        Ma = (i-1)*Ma/i + B/i;
        Ma_inv = (i-1)*Ma_inv/i + inv(B)/i;
        M = sharp(Ma,inv(Ma_inv),0.5);
    elseif strcmp(options.mean,'Inductive')
        M = sharp(M,B,1/i);
    end
    
    %Records the intermediary values
    if options.record
        info.t(i) = toc;
        info.M_rec{i} = M;
    end
end


info.tTot = toc;

end