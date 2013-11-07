clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig3a.m
%
% This program reproduces Fig.3(a) of the paper, which shows MAE vs K
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

K_set     = 2:20;
K_length  = length(K_set);
max_trial = 500;

n = 200;
T = 2;

MAE_LGA  = zeros(max_trial,K_length);
MAE_SBA  = zeros(max_trial,K_length);
MAE_SVD  = zeros(max_trial,K_length);

for i=1:K_length
    fprintf('i = %3g \n', i);
    % Parameter
    K     = K_set(i);
    Q     = K;
    Delta = 0.2;
        
    parfor trial=1:max_trial
        % Construct a Graphon
        w = rand(K,K);
        w = 0.5*(w+w');
        
        % Construct a random graph
        [G1 P_GT1]         = construct_a_graph(w,n,T/2);
        [G2 P_GT2]         = construct_a_graph(w,n/2,T);
        
        clusters_LGA     = estimate_blocks_largest_gap(G1,Q);
        clusters_SBA     = estimate_blocks_directed(G2,Delta);
        [H_LGA P_LGA]    = histogram3D(G1,clusters_LGA);
        [H_SBA P_SBA]    = histogram3D(G2,clusters_SBA);
        P_SVD            = chatterjee_method(G1);
        MAE_LGA(trial,i) = norm(P_LGA(:)-P_GT1(:),1)/numel(P_GT1);
        MAE_SBA(trial,i) = norm(P_SBA(:)-P_GT2(:),1)/numel(P_GT2);
        MAE_SVD(trial,i) = norm(P_SVD(:)-P_GT1(:),1)/numel(P_GT1);
    end
    
end

% load('result_3a');
% save('result_3a');

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

plot(K_set, log10(mean(MAE_SBA)), 'k-o', 'LineWidth', 2); hold on;
plot(K_set, log10(mean(MAE_LGA)), 'k:x', 'LineWidth', 2, 'MarkerSize',8);
plot(K_set, log10(mean(MAE_SVD)), 'k-^', 'LineWidth', 2); hold off;
legend('Proposed', 'Largest Gap', 'USVT','Location','NE');
xlabel('$K$','interpreter','latex');
ylabel('$\log_{10}$(MAE)','interpreter','latex');
grid on;
