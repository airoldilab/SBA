clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig3b.m
%
% This program reproduces Fig.3(b) of the paper, which shows 
% MAE vs Percentage of Missing Links
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

addpath(genpath('../OptSpace_matlab/'));

% Setup Problem
w = [0.8  0.9  0.4  0.5;
     0.1  0.6  0.3  0.2;
     0.3  0.2  0.8  0.3;
     0.4  0.1  0.2  0.9];
n = 200;
T = 2;

p_set     = linspace(0.01,0.2,20);
p_length  = length(p_set);

max_trial = 100;
MAE_MxC  = zeros(max_trial,p_length);
MAE_SBA  = zeros(max_trial,p_length);
MAE_SVD  = zeros(max_trial,p_length);
MAE_LGA  = zeros(max_trial,p_length);

% Main Loop
for i=1:p_length
    fprintf('i = %3g \n', i);
    
    p = p_set(i);

    Delta = 0.2;
    Q = 4;
    
    parfor trial=1:max_trial
        % Observations
        [G1 P_GT1]         = construct_a_graph(w,n,T/2);
        [G2 P_GT2]         = construct_a_graph(w,n/2,T);
        [G3 P_GT3]         = construct_a_graph(w,n/2,T);
        [G4 P_GT4]         = construct_a_graph(w,n,T/2);

        % Masks
        E1 = rand(n,n)>p;
        E2 = rand(n/2,n/2,T)>p;
        E3 = rand(n/2,n/2,T)>p;
        E4 = rand(n,n)>p;
      
        % Chatterjee's Method
        P_SVD            = Method_chatterjee(G1.*E1);
        
        % Our Method
        clusters_SBA     = estimate_blocks_directed(G2.*E2,Delta);
        [~, P_SBA]       = histogram3D_missing(G2.*E2,clusters_SBA,E2);
        
        % Matrix Completion
        P_MxC            = Method_matrix_completion(G3.*E3);
        
        % Chatterjee's Method
        clusters_LGA     = estimate_blocks_largest_gap(G4.*E4,Q);
        [~, P_LGA]       = histogram3D_missing(G4.*E4,clusters_LGA,E4);
        
        MAE_LGA(trial,i) = norm(P_LGA(:)-P_GT4(:),1)/numel(P_GT4);
        MAE_MxC(trial,i) = norm(P_MxC(:)-P_GT3(:),1)/numel(P_GT3);
        MAE_SBA(trial,i) = norm(P_SBA(:)-P_GT2(:),1)/numel(P_GT2);
        MAE_SVD(trial,i) = norm(P_SVD(:)-P_GT1(:),1)/numel(P_GT1);
    end
end

% save('result_fig3b');
% load('result_fig3b');

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

plot(p_set*100, log10(mean(MAE_SBA)), 'k-o', 'LineWidth', 2); hold on;
plot(p_set*100, log10(mean(MAE_LGA)), 'k:x', 'LineWidth', 2, 'MarkerSize',8);
plot(p_set*100, log10(mean(MAE_MxC)), 'k-d', 'LineWidth', 2, 'MarkerSize',8);
plot(p_set*100, log10(mean(MAE_SVD)), 'k-^', 'LineWidth', 2); hold off;
legend('Proposed', 'Largest Gap','OptSpace','USVT', 'Location','SE');
xlabel('$\%$ missing links','interpreter','latex');
ylabel('$\log_{10}$(MAE)','interpreter','latex');
grid on;
