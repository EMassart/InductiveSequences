% Running this code generates Figure 10 of the paper 
% 'Matrix geometric means based on shuffled inductive sequences',
% E. Massart, J. Hendrickx, P.-A. Absil

% Author: E. Massart

clear all;
close all;
clc;
addpath(genpath('methods'));
addpath(genpath('KarcherMean'));

% % -------------------------------------------Computes evolution of the sequences 
n_test = 100;
options1.str = ('para_k10_n3_f1.mat');
validation_parallel(n_test,10,3,1,options1);
options1.str = ('para_k10_n3_f6.mat');
validation_parallel(n_test,10,3,6,options1);

% ---------------------------------------------Plots evolution of the sequences 
linestyle = {'.-', '-s', '-v', '-*'};
names = {'Agg-Ind\_Arith','Agg-Ind\_AH','Agg-Ind\_Ind','Shuff-Ind-Seq'};
files = {'para_k10_n3_f1.mat','para_k10_n3_f6.mat'};
ax = [0 20 10^(-4) 0.2];
N_plot = [10 10]; n_plot = [3 3]; f_plot = [1 6];
col = zeros(4,3);
col(2,:) = [255 102 102];       %rose
col(3,:) = [102 178 255];       %bleu
col(1,:) = [160 160 160];       %gris
col = col/255;
formatSpec = '$N = %d, n = %d, \\kappa  =  10^%d $';
n_meth = 4;

figure;
for i_plot = 1:2
    subplot(1,3,i_plot);
    load(files{i_plot});
    for i_meth = 1:n_meth
        semilogy(1:length(distM(i_meth,:)),distM(i_meth,:),linestyle{i_meth},'Color',col(i_meth,:)); hold on;
    end
    xlabel('$\frac{Number \ of \ \#}{N}$','Interpreter','Latex','Fontsize',12);
    str = sprintf(formatSpec,N_plot(i_plot),n_plot(i_plot),f_plot(i_plot));
    title(str, 'Interpreter', 'Latex');
    axis(ax);
    if i_plot==1
        ylabel('$\mathrm{E_{rel}}(X)$','Interpreter','Latex','Fontsize',12);
    end
end
legend(names{:},'Location','EastOutside');
