clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig2a.m
%
% This program reproduces Fig.2(a) of the paper, which shows MAE vs n
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

addpath(genpath('./OptSpace_matlab/'));

% Graphon
w = [0.8  0.9  0.4  0.5;
     0.1  0.6  0.3  0.2;
     0.3  0.2  0.8  0.3;
     0.4  0.1  0.2  0.9];

% Generate Random Graphs 
n_set    = 100:50:1000; 
n_length = length(n_set);
max_trial = 100;
MAE_LGA  = zeros(max_trial,n_length);
MAE_SBA  = zeros(max_trial,n_length);
MAE_SVD  = zeros(max_trial,n_length);
MAE_MxC  = zeros(max_trial,n_length);

% Main Experiment
for i=1:n_length
    fprintf('i = %3g \n', i);
    n     = n_set(i);
    T     = 2;
    Q     = 4;
    Delta = 0.2;
        
    parfor trial=1:max_trial
        [G1 P_GT1]         = construct_a_graph(w,n,T/2);
        [G2 P_GT2]         = construct_a_graph(w,n/2,T);
        
        clusters_LGA     = estimate_blocks_largest_gap(G1,Q);
        clusters_SBA     = estimate_blocks_directed(G2,Delta);
        [H_LGA P_LGA]    = histogram3D(G1,clusters_LGA);
        [H_SBA P_SBA]    = histogram3D(G2,clusters_SBA);
        P_SVD            = Method_chatterjee(G1);
        P_MxC            = Method_matrix_completion(G2);

        MAE_LGA(trial,i) = norm(P_LGA(:)-P_GT1(:),1)/numel(P_GT1);
        MAE_SBA(trial,i) = norm(P_SBA(:)-P_GT2(:),1)/numel(P_GT2);
        MAE_SVD(trial,i) = norm(P_SVD(:)-P_GT1(:),1)/numel(P_GT1);
        MAE_MxC(trial,i) = norm(P_MxC(:)-P_GT2(:),1)/numel(P_GT2);
    end
end

% save('result_Fig2a');
% load('result_Fig2a');

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

plot(n_set, log10(mean(MAE_SBA)), 'k-o', 'LineWidth', 2); hold on;
plot(n_set, log10(mean(MAE_LGA)), 'k:x', 'LineWidth', 2, 'MarkerSize',8);
plot(n_set, log10(mean(MAE_MxC)), 'k-d', 'LineWidth', 2, 'MarkerSize',8); 
plot(n_set, log10(mean(MAE_SVD)), 'k-^', 'LineWidth', 2); hold off;
legend('Proposed', 'Largest Gap', 'OptSpace','USVT','Location','NE');
xlabel('$n$','interpreter','latex');
ylabel('$\log_{10}$(MAE)','interpreter','latex');
grid on;

