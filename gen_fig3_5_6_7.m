% Running this code generates Figure 3, 5, 6, 7 of the paper 
% 'Matrix geometric means based on shuffled inductive sequences',
% E. Massart, J. Hendrickx, P.-A. Absil

% Author: E. Massart


clear all;
close all;
clc;

addpath(genpath('methods'));
addpath(genpath('KarcherMean'));

% -----------------------------------------------Run the tests
n_test = 30;     %take more datasets when working with few small matrices to remove the (possibly large) variability in computation times arising with such data.
load random_generator.mat;
rng(s);
 
disp('-------------------------------------------------Test 1');
options1.str = ('approx_k3_n3_f1.mat');
validation_opti(n_test,3,3,1,options1,10^(-6));

disp('-------------------------------------------------Test 2');
options1.str = ('approx_k3_n3_f5.mat');
validation_opti(n_test,3,3,5,options1,10^(-6));

n_test = 10;
disp('-------------------------------------------------Test 3');
options1.str = ('approx_k100_n3_f1.mat');
validation_opti(n_test,100,3,1,options1,10^(-6));

disp('-------------------------------------------------Test 4');
options1.str = ('approx_k100_n3_f5.mat');
validation_opti(n_test,100,3,5,options1,10^(-6));

disp('-------------------------------------------------Test 5');
options1.str = ('approx_k30_n100_f1.mat');
validation_opti(n_test,30,100,1,options1,10^(-6));

disp('-------------------------------------------------Test 6');
options1.str = ('approx_k30_n100_f5.mat');
validation_opti(n_test,30,100,5,options1,10^(-6));

disp('-------------------------------------------------Test 7');
options1.str = ('approx_k50_n3_f1.mat');
validation_opti(n_test,50,3,1,options1,10^(-6));

disp('-------------------------------------------------Test 8');
options1.str = ('approx_k10_n100_f1.mat');
validation_opti(n_test,10,100,1,options1,10^(-6));

disp('-------------------------------------------------Test 9');
options1.str = ('approx_k10_n3_f5.mat');
validation_opti(n_test,10,3,5,options1,10^(-6));




%%--------------------------------------------------Plots figure 3
close all;
files = {'approx_k3_n3_f1.mat','approx_k50_n3_f1.mat','approx_k10_n100_f1.mat','approx_k10_n3_f5.mat'};
names = {'Arithmetic','Arithm-Harmo','Log-Euclidean','Cheap','Shuffled Inductive'};
linestyle = {'o','s','v','-^g','.-'};
N_plot = [3 50 10 10]; n_plot = [3 3 100 3]; f_plot = [1 1 1 5];
figure;
col = zeros(5,3);
col(1,:) = [255 102 102];       
col(3,:) = [102 178 255];   
col(2,:) = [255 153 51];     
col(4,:) = [160 160 160];  
col = col/255;
ax = {[0 0.01 10^(-6) 2],[0 0.2 10^(-6) 2],[0 1.2 10^(-6) 2],[0 0.02 10^(-6) 2]};
formatSpec = '$N = %d, n = %d, \\kappa  =  10^%d $';

for i_plot = 1:4
    load(files{i_plot});
    subplot(2,4,i_plot);
    c = 1;
    for i = 1:n_methods
        h(i) = semilogy(qGeo.timeM(c:c+length(qGeo.maxiter{i})-1),qGeo.distM(c:c+length(qGeo.maxiter{i})-1),linestyle{i},'Color',col(i,:)); hold on;
        c = c+length(qGeo.maxiter{i});
    end
    xlabel('Time [s]','Fontsize',9);
    str = sprintf(formatSpec,N_plot(i_plot),n_plot(i_plot),f_plot(i_plot));
    title(str, 'Interpreter', 'Latex');
    axis(ax{i_plot});
    disp('c');
    set(gca, 'YTick', [10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1]);
    if i_plot == 1
        ylabel('$\mathrm{E_{rel}}(X)$','Interpreter','latex','Fontsize',12);
    elseif i_plot == 2 
        legend(h([1 2]),names{[1 2]});
    elseif i_plot == 3
        legend(h([3 4]),names{[3 4]});
    elseif i_plot==4
        legend(h(5),names{5});
    end
end 


%%-----------------------------------------------------Plots figures 5, 6 and 7
clc;
file{1} = {'approx_k3_n3_f1.mat','approx_k3_n3_f5.mat'};
file{2} = {'approx_k100_n3_f1.mat','approx_k100_n3_f5.mat'};
file{3} = {'approx_k30_n100_f1.mat','approx_k30_n100_f5.mat'};
meth_ref = 2;
c = 0.6;
col_1 = [0 133/255 100/255; 1 0 0; 0 0 0; 0 0 1];
n_meth_plot = [9 10:2:18];
col = zeros(length(n_meth_plot),3);
col(2,:) = [0 204 102];             
col(3,:) = [255 153 51];           
col(4,:) = [153 0 153];            
col(5,:) = [51 153 255];            
col(6,:) = [255 102 178];           
col = col/255;
ax{1} = { [0 0.01 10^(-10) 1], [0 0.01 10^(-10) 1], [0 0.01 10^(-10) 1], [0 0.01 10^(-10) 1], [0 0.01 10^(-10) 1], [0 0.01 10^(-10) 1]};
ax{2} = {[0 0.15 10^(-10) 1], [0 0.075 10^(-10) 1], [0 0.075 10^(-10) 1], [0 0.15 10^(-10) 1], [0 0.075 10^(-10) 1], [0 0.075 10^(-10) 1]};
ax{3} = {[0 3 10^(-10) 1],[0 1.5 10^(-10) 1], [0 1.5 10^(-10) 1], [0 3 10^(-10) 1], [0 1.5 10^(-10) 1], [0 1.5 10^(-10) 1]};

for i_plot = 1:3
    figure;
    for i_plot_loc = 1:2   
        load(file{i_plot}{i_plot_loc});
        
        %subplot 1
        subplot(3,3,1+3*(i_plot_loc-1));
        h(1) = semilogy(qGeo.timeM(9:end),qGeo.distM(9:end),'.-','Color',[c c c]); hold on;
        h(2) = semilogy([qGeo.timeM(meth_ref) qGeo.timeM(meth_ref)+reshape(opti.timeM(1,meth_ref,:),1,optionsKarcher.maxiter)],[qGeo.distM(meth_ref) reshape(opti.distM(1,meth_ref,:),1,optionsKarcher.maxiter)],'.-','Color',col_1(1,:)); hold on;
        h(3) = semilogy([qGeo.timeM(9) qGeo.timeM(9)+reshape(opti.timeM(1,9,:),1,optionsKarcher.maxiter)],[qGeo.distM(9) reshape(opti.distM(1,9,:),1,optionsKarcher.maxiter)],'.-','Color',col_1(3,:)); hold on;
        h(4) = semilogy([qGeo.timeM(meth_ref) qGeo.timeM(meth_ref)+reshape(opti.timeM(2,meth_ref,:),1,optionsKarcher.maxiter)],[qGeo.distM(meth_ref) reshape(opti.distM(2,meth_ref,:),1,optionsKarcher.maxiter)],'.-','Color',col_1(4,:)); hold on;
        h(5) = semilogy([qGeo.timeM(9) qGeo.timeM(9)+reshape(opti.timeM(2,9,:),1,optionsKarcher.maxiter)],[qGeo.distM(9) reshape(opti.distM(2,9,:),1,optionsKarcher.maxiter)],'.-','Color',col_1(2,:)); hold on;
        xlabel('Time [s]','Fontsize',9);
        ylabel('$\mathrm{E_{rel}}(X)$','Interpreter','latex','Fontsize',12);
        if i_plot_loc==1
            title('Algorithms comparison','Interpreter','latex','Fontsize',14);
            legend(h(1:3),'Shuffled Inductive seq', 'S.D., X_0 = M_{AH}', 'S.D., X_0 = X^S_{N}');
        else
            legend(h(4:5),'LRBFGS, X_0 = M_{AH}','LRBFGS, X_0 = X^S_{N}');
        end
        set(gca,'ytick',[10^(-8) 10^(-6) 10^(-4) 10^(-2) 1]);
        axis(ax{i_plot}{1+3*(i_plot_loc-1)});
        
        %subplot 2
        subplot(3,3,2+3*(i_plot_loc-1));
        names = cell(1,length(n_meth_plot));
        linestyle = cell(1,length(n_meth_plot));
        names{1} = strcat('k = ',num2str(1));
        for i_meth = 1:length(n_meth_plot)-1
            names{i_meth+1} = strcat('k = ',num2str(2*i_meth));
            linestyle{i_meth} = '.-';           
        end
        semilogy(qGeo.timeM(9:end),qGeo.distM(9:end),'.-','Color',[c c c]); hold on;
        for i = 1:length(n_meth_plot)
            h(i) = semilogy([qGeo.timeM(n_meth_plot(i)) qGeo.timeM(n_meth_plot(i))+reshape(opti.timeM(1,n_meth_plot(i),:),1,optionsKarcher.maxiter)],[qGeo.distM(n_meth_plot(i)) reshape(opti.distM(1,n_meth_plot(i),:),1,optionsKarcher.maxiter)],'.-','Color',col(i,:)); hold on;
        end
        xlabel('Time [s]','Fontsize',9);
        if i_plot_loc==1
            title('Initialization SD','Interpreter','Latex','Fontsize',14);
        else
            legend('Shuffled_ Inductive seq',names{1:2});
        end
        set(gca,'ytick',[10^(-8) 10^(-6) 10^(-4) 10^(-2) 1]);
        axis(ax{i_plot}{2+3*(i_plot_loc-1)});
        
        %subplot 3
        subplot(3,3,3+3*(i_plot_loc-1));
        semilogy(qGeo.timeM(9:end),qGeo.distM(9:end),'.-','Color',[c c c]); hold on;
        for i = 1:length(n_meth_plot)
            h(i) = semilogy([qGeo.timeM(n_meth_plot(i)) qGeo.timeM(n_meth_plot(i))+reshape(opti.timeM(2,n_meth_plot(i),:),1,optionsKarcher.maxiter)],[qGeo.distM(n_meth_plot(i)) reshape(opti.distM(2,n_meth_plot(i),:),1,optionsKarcher.maxiter)],'.-','Color',col(i,:)); hold on;
        end
        xlabel('Time [s]','Fontsize',9);
        if i_plot_loc==1
            title('Initialization LRBFGS','Interpreter','Latex','Fontsize',14);
        else
            legend(h([3:6]),names{3:6});
        end
        set(gca,'ytick',[10^(-8) 10^(-6) 10^(-4) 10^(-2) 1]);
        axis(ax{i_plot}{3+3*(i_plot_loc-1)});
    end
end


