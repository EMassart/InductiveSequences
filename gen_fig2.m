% Running this code generates Figure 2 of the paper 
% 'Matrix geometric means based on shuffled inductive sequences',
% E. Massart, J. Hendrickx, P.-A. Absil

% Author: E. Massart
% 
clear all;
close all;
clc;

addpath(genpath('methods'));
addpath(genpath('KarcherMean'));

% ---------------------------------------First row : computes the evolution of the sequences 
n_test = 1000;
options1.maxiter = 10;
options1.str = ('shuffle_k3_n3_f1.mat');
graphe_comparison_shuffle(n_test,3,3,1,options1);
options1.str = ('shuffle_k50_n3_f1.mat');
graphe_comparison_shuffle(n_test,50,3,1,options1);
options1.str = ('shuffle_k10_n100_f1.mat');
graphe_comparison_shuffle(n_test,10,100,1,options1);
options1.str = ('shuffle_k10_n3_f6.mat');
graphe_comparison_shuffle(n_test,10,3,6,options1);

% ---------------------------------------First row: plots the graphs
n_meth = 4;
names = {'Circular','Random without replacement','Random with replacement','Shuffled'};
files = {'shuffle_k3_n3_f1.mat','shuffle_k50_n3_f1.mat','shuffle_k10_n100_f1.mat','shuffle_k10_n3_f6.mat'};


ax = [0 10 5*10^(-4) 1];
N_plot = [3 50 10 10]; n_plot = [3 3 100 3]; f_plot = [1 1 1 6];
linestyle = {'-v','-s','-^','-o'};
linestyleSTD = {'-.','-.','-.','-.'};

col = zeros(4,3);
col(2,:) = [255 102 102];       
col(3,:) = [102 178 255];       
col(1,:) = [160 160 160];       
col = col/255;
formatSpec = '$N = %d, n = %d, \\kappa  =  10^%d $';

figure;
for i = 1:4
    subplot(2,4,i);
    load(files{i});
    STD = reshape(std(dist,1),n_methods,options.maxiter);
    for i_meth = 1:n_meth
        semilogy(1:options1.maxiter,distM(i_meth,:)+STD(i_meth,:),linestyleSTD{i_meth},'Color',col(i_meth,:)); hold on;
        semilogy(1:options1.maxiter,distM(i_meth,:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
        semilogy(1:options1.maxiter,distM(i_meth,:)-STD(i_meth,:),linestyleSTD{i_meth},'Color',col(i_meth,:)); hold on;
    end
    xlabel('$k$','Interpreter','Latex','Fontsize',12);
    str = sprintf(formatSpec,N_plot(i),n_plot(i),f_plot(i));
    title(str, 'Interpreter', 'Latex');
    axis(ax);
    if i>1
        set(gca, 'YTickLabel', []);
    else
        ylabel('$\mathrm{E_{rel}}(X_{kN})$','Interpreter','Latex','Fontsize',12);
    end
end


% ----------------------------------------Second row: computes the evolution with params
clear all; clc;
n_test = 100;
options.maxiter = 11;
n_meth = 4;

param{1} = {[3 5 10:10:50], 3, 1};
param{2} = {10, [3 10:10:100], 1};
param{3} = {10, 3, [1:6]};

for i = 1:3
    fprintf('-----------------------------------------------------Test number: %d \n',i);
    dist{i} = zeros(n_meth,length(param{i}{i}));
    distSTD{i} = zeros(n_meth,length(param{i}{i}));
    p = zeros(1,3);
    p([1:i-1, i+1:end]) = [param{i}{[1:i-1, i+1:end]}];
    for i_loc = 1:length(param{i}{i})
        fprintf('-----------------------------------------------------Parameter: %d \n',param{i}{i}(i_loc));
        p(i) = param{i}{i}(i_loc);
        p = p
        [~,dist{i}(:,i_loc),distSTD{i}(:,i_loc)] = graphe_comparison_shuffle(n_test,p(1),p(2),p(3),options);
    end
end
save('Comparison_Shuffle_param.mat');


%% ---------------------------------------Second row : plots the evolution with params
load('Comparison_Shuffle_param.mat');
names = {'Circular','Rand. without repl.','Rand. with repl.','Shuffled'};
linestyle = {'-v','-s','-^','-o'};

col = zeros(4,3);
col(2,:) = [255 102 102];       
col(3,:) = [102 178 255];       
col(1,:) = [160 160 160];      
col = col/255;
ax = {[0 50 5*10^(-4) 1],[0 100 5*10^(-4) 1],[10 10^6 5*10^(-4) 1]};
str = { '$n = 3, \\kappa = 10 $', '$N = 10, \\kappa = 10 $', '$N = 10, n = 3$'};
str_label = {'$N$','$n$','$\\kappa$'};
for i = 1:2;
    subplot(2,4,i+4);
    for i_meth = 1:n_meth
        semilogy(param{i}{i},dist{i}(i_meth,:)+distSTD{i}(i_meth,:),'-.','Color',col(i_meth,:)); hold on;
        semilogy(param{i}{i},dist{i}(i_meth,:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
        semilogy(param{i}{i},dist{i}(i_meth,:)-distSTD{i}(i_meth,:),'-.','Color',col(i_meth,:)); hold on;
    end
    xlabel(sprintf(str_label{i}),'Interpreter','Latex','Fontsize',12);
    ylabel('$\mathrm{E_{rel}}(X_{10N})$','Interpreter','Latex','Fontsize',12);
    axis(ax{i});
    title(sprintf(str{i}), 'Interpreter', 'Latex');
end
set(gca, 'YTickLabel', []);
subplot(2,4,7);
for i_meth = 1:n_meth
    loglog(10.^(param{3}{3}),dist{3}(i_meth,:)+distSTD{3}(i_meth,:),'-.','Color',col(i_meth,:)); hold on;
    h(i_meth) = loglog(10.^(param{3}{3}),dist{3}(i_meth,:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
    loglog(10.^(param{3}{3}),dist{3}(i_meth,:)-distSTD{3}(i_meth,:),'-.','Color',col(i_meth,:)); hold on;
end
xlabel('$\kappa$','Interpreter','Latex','Fontsize',12);
axis(ax{3});
title(str{3}, 'Interpreter', 'Latex');
set(gca, 'YTickLabel', []);
set(gca, 'XTick', [ 10^2 10^4 10^6 ]);
legend(h(:),names{:},'Location','EastOutside');
