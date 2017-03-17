% Running this code generates Figure 4, 8, 9 of the paper 
% 'Matrix geometric means based on shuffled inductive sequences',
% E. Massart, J. Hendrickx, P.-A. Absil

% Author: E. Massart


%% --------------------------------------------------Test to generate figures 4, 8 and 9

clear all; close all; clc;
addpath(genpath('methods'));
addpath(genpath('KarcherMean'));

n_test = 50;                 
n_meth = 18;
load random_generator.mat;
rng(s);

param{1} = {[3 5 10:10:50], 3, 1};
param{2} = {10, [3 10:10:100], 1};
param{3} = {10, 3, [1:6]};

for i = 1:3
    fprintf('-----------------------------------------------------Test number: %d \n',i);
    dist{i} = zeros(n_meth,length(param{i}{i}));
    time{i} = zeros(n_meth,length(param{i}{i}));
    nIter{i} = zeros(n_meth,length(param{i}{i}));
    nIter2{i} = zeros(n_meth,length(param{i}{i}));
    tStop{i} = zeros(n_meth,length(param{i}{i}));
    tStop2{i} = zeros(n_meth,length(param{i}{i}));
    p = zeros(1,3);
    p([1:i-1, i+1:end]) = [param{i}{[1:i-1, i+1:end]}];
    for i_loc = 1:length(param{i}{i})
        fprintf('-----------------------------------------------------Parameter: %d \n',param{i}{i}(i_loc));
        p(i) = param{i}{i}(i_loc);
        p = p
        [qGeo,opti] = validation_opti(n_test,p(1),p(2),p(3),struct(),10^(-6));
        dist{i}(:,i_loc) = qGeo.distM(1:18);
        time{i}(:,i_loc) = qGeo.timeM(1:18);
        nIter{i}(:,i_loc) = opti.kStopM(1,:)';
        nIter2{i}(:,i_loc) = opti.kStopM(2,:)';
        tStop{i}(:,i_loc) = opti.tStopM(1,:)';
        tStop2{i}(:,i_loc) = opti.tStopM(2,:)';
    end
    save(strcat('A',num2str(i),'.mat'));
end



%% ---------------------------------------------------Plots evolution params
names = {'Arithmetic','Arithm-Harmo','Log-Euclidean','Cheap: k_{Ch} = 1','Shuff. Inductive: k = 1 (= M_{Ind})','Shuff. Inductive: k = 10'};
linestyle = {'-o','-s','-v','-^','-o','-*'};
n_meth_plot = [1:4 9 12];
col = zeros(6,3);
col(1,:) = [255 102 102];      
col(3,:) = [102 178 255];      
col(2,:) = [255 153 51];       
col(4,:) = [160 160 160];     
col = col/255;
ax1 = {[0 50 10^(-3) 2],[0 100 10^(-3) 2],[10 10^6 10^(-3) 2]};
ax2 = {[0 50 10^(-5) 0.2],[0 100 10^(-5) 0.2],[10 10^6 10^(-5) 0.2]};
str = { '$n = 3, \\kappa = 10 $', '$N = 10, \\kappa = 10 $', '$N = 10, n = 3$'};
str_label = {'$N$','$n$','$\\kappa$'};
for i = 1:2;
    load(strcat('A',num2str(i),'.mat'));
    subplot(2,4,i);
    for i_meth = 1:length(n_meth_plot)
        h(i_meth) = semilogy(param{i}{i},dist{i}(n_meth_plot(i_meth),:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
    end
    xlabel(sprintf(str_label{i}),'Interpreter','Latex','Fontsize',12);
    ylabel('$\mathrm{E_{rel}}(X)$','Interpreter','Latex','Fontsize',12);
    axis(ax1{i});
    title(sprintf(str{i}), 'Interpreter', 'Latex');
    set(gca,'ytick',[10^(-3) 10^(-2) 10^(-1)  1])
    
    subplot(2,4,i+4);
    for i_meth = 1:length(n_meth_plot)
        h(i_meth) = semilogy(param{i}{i},time{i}(n_meth_plot(i_meth),:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
    end
    xlabel(sprintf(str_label{i}),'Interpreter','Latex','Fontsize',12);
    ylabel('Time[s]','Interpreter','Latex','Fontsize',12);
    axis(ax2{i});
    set(gca,'ytick',[10^(-5) 10^(-3)  10^(-1)])
end

subplot(2,4,3);
load('A3.mat');
for i_meth = 1:length(n_meth_plot)
    loglog(10.^(param{3}{3}),dist{3}(n_meth_plot(i_meth),:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
end
xlabel(sprintf(str_label{3}),'Interpreter','Latex','Fontsize',12);
axis(ax1{3});
title(sprintf(str{i}), 'Interpreter', 'Latex');
set(gca,'ytick',[10^(-3) 10^(-2) 10^(-1)  1])
set(gca, 'XTick', [10 10^3 10^5]);
% 
subplot(2,4,7);
for i_meth = 1:length(n_meth_plot)
     loglog(10.^(param{3}{3}),time{3}(n_meth_plot(i_meth),:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
end
xlabel(sprintf(str_label{3}),'Interpreter','Latex','Fontsize',12);
legend(h(1:length(n_meth_plot)),names{:},'Location','EastOutside');
axis(ax2{3});
set(gca,'ytick',[10^(-5) 10^(-3) 10^(-1)])
set(gca, 'XTick', [10 10^3 10^5]);


%% -----------------------------------------------Plots SD initialization
% to generate figure 9 (LRBFGS initialization), simply replace tStop by tStop2

clear all;
close all;
clc;

files = {'A1.mat', 'A2.mat', 'A3.mat'};
figure; 
n_meth_plot = [1:4 9 12];
col = zeros(6,3);
col(1,:) = [255 102 102];       
col(3,:) = [102 178 255];    
col(2,:) = [255 153 51];   
col(4,:) = [160 160 160];      
col = col/255;
linestyle = {'-o','-s','-v','-^','-o','-*'};
names = {'Arithmetic','Arithm-Harmo','Log-Euclidean','Cheap: k = 1','Shuffled Inductive: k = 1','Shuffled Inductive: k = 4'};
ax = {[0 50 10^(-3) 0.3], [0 100 10^(-3) 0.3],[10 10^6 10^(-3) 0.3]};

str = { '$n = 3, \\kappa = 10 $', '$N = 10, \\kappa = 10 $', '$N = 10, n = 3$'};
str_label = {'$N$','$n$','$\\kappa$'};

for i_plot = 1:2
    load(files{i_plot});
    subplot(1,4,i_plot);
    for i = 1:length(n_meth_plot)
        h(i) = semilogy(param{i_plot}{i_plot},tStop{i_plot}(n_meth_plot(i),:)+time{i_plot}(n_meth_plot(i),:),linestyle{i},'Color',col(i,:)); hold on;
    end
    xlabel(sprintf(str_label{i_plot}),'Interpreter','Latex','Fontsize',12);
    if i_plot==1
        ylabel('time [s]','Fontsize',12);
    end
    title(sprintf(str{i_plot}), 'Interpreter', 'Latex');
    axis(ax{i_plot});
end

load(files{3});
subplot(1,4,3);
for i_meth = 1:length(n_meth_plot)
    loglog(10.^(param{3}{3}),tStop{3}(n_meth_plot(i_meth),:)+time{3}(n_meth_plot(i_meth),:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
end
xlabel(sprintf(str_label{3}),'Interpreter','Latex','Fontsize',12);
axis(ax{3});
legend(names{:},'Location','EastOutside');
title(sprintf(str{i}), 'Interpreter', 'Latex');
set(gca, 'XTick', [10 10^3 10^5]);

