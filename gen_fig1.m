% Running this code generates Figure 1 of the paper 
% 'Matrix geometric means based on shuffled inductive sequences',
% E. Massart, J. Hendrickx, P.-A. Absil

% Author: E. Massart

clear all;
close all;
clc;

addpath(genpath('methods'));
addpath(genpath('KarcherMean'));

N = 3;
n = 10;
n_test = 1000;
k = 10;
optionsKarcher.algo = 'sum';
optionsKarcher.maxiter = 1000;
optionsKarcher.stop = 1;
dist_rec = zeros(n_test,N,k*N);

for i_test = 1:n_test
    
    if mod(i_test,100)==0 fprintf('i_test = %d \n',i_test); end
    %data
    problem.size = n;
    problem.number = N;
    [ A ] = gen_mat( problem );
    
    %Karcher exact
    meanKarcher = karcher(A, optionsKarcher);
    
    d = zeros(1,N);
    for i_mat = 1:N
        d(i_mat) = dist_mat(A{i_mat},meanKarcher);
    end
    
    X = A{1};
    alpha = 1;
    for i_t = 1:k*N
        i_mat = mod(i_t,N);
        if (i_mat == 0), i_mat = N; end
        X = sharp(X,A{i_mat},1/alpha);
        alpha = alpha + 1;
        for i_mat = 1:N
            dist_rec(i_test,i_mat,i_t) = dist_mat(A{i_mat},X)/d(i_mat);
        end
    end
end
distM = mean(dist_rec,1);


%Plots
c = 0.5;
k = 6;
distM = distM(1,:,1:N*k);
distM_bis = distM(1,:,N:N:k*N);
figure;
hold on;
plot(1:k*N,reshape(distM(1,1,:),1,N*k),'o','Color',[c c c]);
plot(1:k*N,reshape(distM(1,2,:),1,N*k),'x','Color',[c c c]);
plot(1:k*N,reshape(distM(1,3,:),1,N*k),'v','Color',[c c c]);
plot(N:N:k*N,reshape(distM_bis(1,1,:),1,k),'-ok');
plot(N:N:k*N,reshape(distM_bis(1,2,:),1,k),'-xk');
plot(N:N:k*N,reshape(distM_bis(1,3,:),1,k),'-vk');
set(gca,'XTick',0:N:k*N);
xlabel('$k$','Interpreter','Latex','Fontsize',15);
ylabel('$\frac{\delta(X_k,A_i)}{\delta(K,A_i)}$','Interpreter','Latex','Fontsize',15);
legend('i = 1','i = 2','i = 3');
save('figure1.mat');


