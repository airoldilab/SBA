clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig4.m
%
% This program reproduces Fig.4 of the paper, which shows MAE vs n for two
% types of graphons
% 
% Case 1: w(u,v) = 1/(1 + exp(-50(u^2 + v^2)))
% Case 2: w(u,v) = uv
% 
% The proposed algorithm is compared with
% (1) Universal Single Value Thresholding (by S. Chatterjee, 2012)
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

% Random Graphs 
n_set    = 100:50:1000; 
n_length = length(n_set);
max_trial = 100;
MAE_SBA  = zeros(max_trial,n_length);
MAE_SVD  = zeros(max_trial,n_length);

caseid = 1; % Switch between case 1 and case 2

for ii=1:n_length
    fprintf('ii = %3g \n', ii);
    n     = n_set(ii);
    T     = 2;
    Delta = 0.1;
    
    parfor trial=1:max_trial
        
        I = linspace(0,1,n);
        J = linspace(0,1,n);
        P_GT = zeros(n,n);
        for i=1:n
            for j=1:n
                if caseid==1
                    P_GT(i,j) = 1/(1+exp(-50*(I(i)+J(j))));
                end
                
                if caseid==2
                    P_GT(i,j) = I(i)*J(j);
                end
            end
        end
        P_GT1 = P_GT;
        P_GT2 = P_GT(1:2:n,1:2:n);
        
        G1 = construct_a_graph_from_P(P_GT1,n,T/2);
        G2 = construct_a_graph_from_P(P_GT2,n/2,T);
        
        clusters_SBA     = estimate_blocks_directed(G2,Delta);
        [H_SBA P_SBA]    = histogram3D(G2,clusters_SBA);
        P_SVD            = chatterjee_method(G1);
        
        MAE_SBA(trial,ii) = norm(P_SBA(:)-P_GT2(:),1)/numel(P_GT2);
        MAE_SVD(trial,ii) = norm(P_SVD(:)-P_GT1(:),1)/numel(P_GT1);
    end
end



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
plot(n_set, log10(mean(MAE_SVD)), 'k-^', 'LineWidth', 2); hold off;
legend('Proposed', 'USVT','Location','NE');
xlabel('$n$','interpreter','latex');
ylabel('$\log_{10}$(MAE)','interpreter','latex');
grid on;



