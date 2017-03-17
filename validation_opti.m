function [qGeo,opti] = validation_opti(n_test,k,n,f,options_glob,tol)
% Generates the validation plots
% Author: E. Massart

methods = {'arithm','arithm_harmo','explog','cheap','Inductive'};
names = {'Arithmetic','Arithm-Harmo','Log-Euclidean','Cheap','Shuffled Inductive'};

if ~exist('options_glob','var')
    options_glob = struct();
end

n_methods = length(methods);

%params for the first test: computes quasi_geometric means
qGeo = struct();
qGeo.maxiter = cell(1,length(methods));
qGeo.maxiter(1:3) = {1};
qGeo.maxiter(4) = {1:5};
qGeo.maxiter(5) = {1:10};
qGeo.n_methods = length(cat(2,qGeo.maxiter{:}));         % total number of quasi-geometric means
qGeo.time = zeros(n_test, qGeo.n_methods);
qGeo.dist = zeros(n_test, qGeo.n_methods);

%params for the second test: run opti methods using qGeo means as initial point
opti = struct();
opti.maxiter = cell(1,length(methods));
opti.maxiter(1:3) = {1};
opti.maxiter(4) = {1:5};
opti.maxiter(5) = {1:10};
opti.n_methods = length(cat(2,opti.maxiter{:}));

optionsKarcher.maxiter = 500;
optionsLRBFGS.maxiter = optionsKarcher.maxiter;         %choose those two parameters to be equal to allow to record the results in a same array.
optionsKarcher.stop = 1;

optionsLRBFGS.stop = 5;
optionsLRBFGS.err = 10^(-12);
opti.dist = zeros(2, n_test, opti.n_methods, optionsKarcher.maxiter);
opti.time = zeros(2, n_test, opti.n_methods, optionsKarcher.maxiter);
opti.kStop = zeros(2, n_test, opti.n_methods);
opti.tStop = zeros(2, n_test, opti.n_methods);

d_ref = zeros(1,n_test);

problem.size = n;
problem.number = k;
problem.cond = f;

optionsKarcher_ref.maxiter = 1000;
optionsKarcher_ref.stop = 1;
i_test = 1;
while i_test <= n_test
    
    % problem generation
    if mod(i_test,10)==0 fprintf('test number %d \n',i_test); end

    [ A ] = gen_mat( problem );
    Kmean = karcher(A,optionsKarcher_ref);
    
    d_loc = zeros(1,k);
    for i_mat = 1:k
        d_loc(i_mat) = dist_mat(A{i_mat},Kmean);                
    end
    d_ref(i_test) = mean(d_loc);
    
    M_start = cell(1,opti.n_methods);
    qGeo.c = 1;                                         %counter to fill qGeo.dist and qGeo.time (all quasi geometric means)
    opti.c = 1;                                         %counter to fill M_start (only quasi geom means used as initial point for opti methods)
    
    % first test: computes quasi geometric means
    for i_meth = 1:n_methods
        options.maxiter = max(qGeo.maxiter{i_meth});
        options.record = 1;
        [M,info] = eval([methods{i_meth},'( A, options );']);
        if isfield(info,'M_rec')
            for i = 1:length(qGeo.maxiter{i_meth})
                if ~isempty(find(opti.maxiter{i_meth}==qGeo.maxiter{i_meth}(i),1))
                    M_start{opti.c} = info.M_rec{qGeo.maxiter{i_meth}(i)};
                    opti.c = opti.c+1;
                end
                qGeo.dist(i_test,qGeo.c) = dist_mat(info.M_rec{qGeo.maxiter{i_meth}(i)},Kmean)./d_ref(i_test);
                qGeo.time(i_test,qGeo.c) = info.t(qGeo.maxiter{i_meth}(i));
                qGeo.c = qGeo.c+1;
            end
        else
            M_start{opti.c} = M;
            opti.c = opti.c+1;
            qGeo.dist(i_test,qGeo.c) = dist_mat(M,Kmean)./d_ref(i_test);
            qGeo.time(i_test,qGeo.c) = info.tTot;
            qGeo.c = qGeo.c+1;
        end
    end
    
    
    % Test 2 : SD method
    stop_SD = 0;
    optionsKarcher.MStart =  eye(problem.size);
    karcher( A, optionsKarcher );  
    for i_meth = 1:opti.n_methods
        optionsKarcher.MStart = M_start{i_meth} ;
        [~,info] = karcher( A, optionsKarcher );
        for i = 1:length(info.M_rec)
            opti.dist(1,i_test,i_meth,i) = dist_mat(info.M_rec{i},Kmean)./d_ref(i_test);
            opti.time(1,i_test,i_meth,i) = info.t(i);
        end
        if length(info.M_rec)<optionsKarcher.maxiter                            %if the method converges before reaching the max number of iterations, the remaining zeros are replaced by the last accuracy obtained.
            opti.dist(1,i_test,i_meth,length(info.M_rec)+1:end) = opti.dist(1,i_test,i_meth,length(info.M_rec));
            opti.time(1,i_test,i_meth,length(info.M_rec)+1:end) = opti.time(1,i_test,i_meth,length(info.M_rec));
        end
        if ~isempty(find(reshape(opti.dist(1,i_test,i_meth,:),1,optionsKarcher.maxiter)<tol*d_ref(i_test),1))
            opti.kStop(1,i_test,i_meth) = find(reshape(opti.dist(1,i_test,i_meth,:),1,optionsKarcher.maxiter)<tol*d_ref(i_test),1);
            opti.tStop(1,i_test,i_meth) = opti.time(1,i_test,i_meth,opti.kStop(1,i_test,i_meth));
        else 
            i_test = i_test-1;
            stop_SD = 1;
            fprintf('Convergence SD failed: max number of iteration reached before reaching tol of %14.8e, method %d \n',tol,i_meth);
            break;
        end
    end
    
    if ~stop_SD                                                             %if SD converges properly on the current dataset, try to run the LRBFGS method using the same data points. Otherwise, generate a new data set.
        % Test 2 : LRBFGS method
        optionsLRBFGS.MStart = eye(problem.size);
        LRBFGS_toolbox( A, optionsLRBFGS );
        for i_meth = 1:opti.n_methods
            optionsLRBFGS.MStart = M_start{i_meth};
            [~,info] = LRBFGS_toolbox( A, optionsLRBFGS );
            if info.status >= 0
                for i = 2:length(info.M_rec)
                    opti.dist(2,i_test,i_meth,i-1) = dist_mat(info.M_rec{i}.U,Kmean)./d_ref(i_test);
                    opti.time(2,i_test,i_meth,i-1) = info.t(i);
                end
                if length(info.M_rec)-1<optionsLRBFGS.maxiter                   %if the method converges before reaching the max number of iterations, the remaining zeros are replaced by the last accuracy obtained.
                    opti.dist(2,i_test,i_meth,length(info.M_rec):end) = opti.dist(2,i_test,i_meth,length(info.M_rec)-1);
                    opti.time(2,i_test,i_meth,length(info.M_rec):end) = opti.time(2,i_test,i_meth,length(info.M_rec)-1);
                end
                if ~isempty(find(reshape(opti.dist(2,i_test,i_meth,:),1,optionsLRBFGS.maxiter)<tol*d_ref(i_test),1))
                    opti.kStop(2,i_test,i_meth) = find(reshape(opti.dist(2,i_test,i_meth,:),1,optionsLRBFGS.maxiter)<tol*d_ref(i_test),1);
                    opti.tStop(2,i_test,i_meth) = opti.time(2,i_test,i_meth,opti.kStop(2,i_test,i_meth));
                else
                    i_test = i_test-1;
                    fprintf('Convergence LRBFGS failed: max number of iteration reached before reaching tol of %14.8e, method %d \n',tol,i_meth);
                    break;
                end
            else
                i_test = i_test-1;
                fprintf('Convergence LRBFGS failed, method %d \n',i_meth);
                break;
            end
        end
    end
    i_test = i_test+1;
end

qGeo.timeM = reshape(mean(qGeo.time,1),1,qGeo.n_methods);
qGeo.distM = reshape(mean(qGeo.dist,1),1,qGeo.n_methods);
opti.distM = reshape(mean(opti.dist,2),2,opti.n_methods,optionsKarcher.maxiter);
opti.timeM = reshape(mean(opti.time,2),2,opti.n_methods,optionsKarcher.maxiter);
opti.kStopM =  reshape(mean(opti.kStop,2),2,opti.n_methods);
opti.tStopM =  reshape(mean(opti.tStop,2),2,opti.n_methods);

if isfield(options_glob,'str')
    save(options_glob.str);
end

end
