function [tTot,dTot,dSTDTot] =  graphe_comparison_shuffle(n_test,k,n,f,options)
% Author: E. Massart

if ~exist('options','var')
    options = struct();
end

if ~isfield(options,'maxiter')
    options.maxiter = 100;
end

methods = {'no','random','random_replacement','shuffling'};
n_methods = length(methods);
time = zeros(n_test, n_methods, options.maxiter);
dist = zeros(n_test, n_methods, options.maxiter);
d_ref = zeros(1,n_test);

for i_test = 1:n_test
    
    if mod(i_test,10)==0 fprintf('Set number %d \n',i_test); end
    
    %generate data
    %[Kmean, A, ~, ~] = constructed_data(k,n,f);
    problem.size = n;
    problem.number = k;
    problem.cond = f;
    [ A ] = gen_mat( problem );
    optionsKarcher_ref.maxiter = 1000;
    optionsKarcher_ref.stop = 1;
    Kmean = karcher(A,optionsKarcher_ref);
    
    d_loc = zeros(1,k);
    for i_mat = 1:k
        d_loc(i_mat) = dist_mat(A{i_mat},Kmean);
    end
    d_ref(i_test) = mean(d_loc);
    
    optionsMeth.record = 1;
    optionsMeth.maxiter = options.maxiter;
    for i_meth = 1:n_methods
        optionsMeth.shuffling = methods{i_meth};
        [~,info] = Inductive( A, optionsMeth );
        for i_iter = 1:optionsMeth.maxiter
            time(i_test,i_meth,i_iter) = info.t(i_iter);
            dist(i_test,i_meth,i_iter) = dist_mat(Kmean,info.M_rec{i_iter})./d_ref(i_test);
        end
    end
end

timeM = reshape(mean(time,1),n_methods,options.maxiter);
distM = reshape(mean(dist,1),n_methods,options.maxiter);
distSTD = reshape(std(dist,1),n_methods,options.maxiter);
tTot = timeM(:,options.maxiter);
dTot = distM(:,options.maxiter);
dSTDTot = distSTD(:,options.maxiter);

if isfield(options,'str')
    save(options.str);
end

end