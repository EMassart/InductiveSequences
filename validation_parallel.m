function [distM] = validation_parallel(n_test,k,n,f,options_glob)
% Generates the validation plots
% Author: E. Massart


methods = {'Arith','AH','Inductive'};

if ~exist('options_glob','var')
    options_glob = struct();
end

n_methods = length(methods)+1;
options.maxiter = 20;
dist = zeros(n_test, n_methods, options.maxiter);

for i_test = 1:n_test
    
    if mod(i_test,10)==0 fprintf('test number %d \n',i_test); end
    
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
    
    for i_meth =  1:n_methods-1
        options.record = 1;
        options.mean = methods{i_meth};
        options.shuffling = 'shuffling';
        [~,info] = Inductive_parallel( A, options );
        for i = 1:options.maxiter
            dist(i_test,i_meth,i) = dist_mat(info.M_rec{i},Kmean)./d_ref(i_test);
        end
    end
    i_meth = i_meth + 1;
    options2.record = 1;
    options2.maxiter = options.maxiter;
    options2.shuffling = 'shuffling';
    [~,info] = Inductive( A, options2 );
    for i = 1:options2.maxiter
        dist(i_test,i_meth,i) = dist_mat(info.M_rec{i},Kmean)./d_ref(i_test);
    end
    
end
distM = reshape(mean(dist,1),n_methods,options.maxiter);

if isfield(options_glob,'str')
    save(options_glob.str);
end

end
