clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig2b.m
%
% This program reproduces Fig.2(b) of the paper, which shows MAE vs T
% 
% The proposed algorithm is compared with
% (1) Large Gap Method (by A. Channarond, J. Daudin, and S. Robin, 2012)
% (2) Universal Single Value Thresholding (by S. Chatterjee, 2012)
% (3) Matrix Completion (by R.H. Keshavan, A.Montanari, and S. Oh, 2010)
%
%
% Reference
% E. M. Airoldi, T. B. Costa, S. H. Chan, "Stochastic blockmodel approximation of a graphon:
% Theory and consistent estimation", Advances in Neural Information
% Processing Systems, 2013
%
% 
% copy-right 2013
% Harvard University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Graphon
w = [0.8  0.9  0.4  0.5;
     0.1  0.6  0.3  0.2;
     0.3  0.2  0.8  0.3;
     0.4  0.1  0.2  0.9];

% Random Graphs 
T_set    = 2:2:40;
T_length = length(T_set);
max_trial = 100;
MAE_SBA  = zeros(max_trial,T_length);

for i=1:T_length
    fprintf('i = %3g \n', i);
    T     = T_set(i);
    n     = 200;
    Q     = 3;
    Delta = 0.2;
        
    parfor trial=1:max_trial
        [G2 P_GT2]       = construct_a_graph(w,n,T);
        clusters_SBA     = estimate_blocks_directed(G2,Delta);
        [H_SBA P_SBA]    = histogram3D(G2,clusters_SBA);
        MAE_SBA(trial,i) = norm(P_SBA(:)-P_GT2(:),1)/numel(P_GT2);
    end
end

% save('result_fig2b');
% load('result_fig2b');

figure(1);
fontsize = 12;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
fontname = 'Times New Roman';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
fontweight = 'normal';
set(0,'defaultaxesfontweight',fontweight);
set(0,'defaulttextfontweight',fontweight);

plot(T_set, log10(mean(MAE_SBA)), 'k-o', 'LineWidth', 2);
legend('Proposed','Location','NE');
xlabel('$2T$','interpreter','latex');
ylabel('$\log_{10}$(MAE)','interpreter','latex');
grid on;

