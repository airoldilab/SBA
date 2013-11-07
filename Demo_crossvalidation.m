clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo file for executing the cross validation
%
% Stanley Chan @ Harvard
% Apr 23, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = [0.7 0.2; 0.2 0.5];
n = 200;
T = 2;

[G P_GT] = construct_a_graph(w,n,T);

Delta_Set = linspace(0.01,0.4,20);
max_trial = 5;
mse       = zeros(length(Delta_Set),max_trial);
J         = zeros(length(Delta_Set),max_trial);
Bins      = zeros(length(Delta_Set),max_trial);

parfor i=1:length(Delta_Set);
    fprintf('i = %3g \n', i);
    
    for trial = 1:max_trial
        Delta = Delta_Set(i);
        
        % Estimate the blocks
        B = estimate_blocks_directed(G,Delta);
        
        % Calculate MSE
        [H, P_est]   = histogram3D(G,B);
        mse(i,trial) = compute_mse(P_GT,P_est);
        
        % Calculate Cross Validation
        m = length(B);
        p = zeros(m,1);
        for j=1:m
            p(j) = length(B{j})/n;
        end
        h = 1/m;
        J(i, trial)  = 2/(h*(n-1)) - (n+1)/(h*(n-1))*sum(p.^2);
        
        % Calculate Num of Bins
        Bins(i, trial) = m;
    end
end
[val pos] = min(mean(J,2));

P_char    = Method_chatterjee(G);
rmse_char = sqrt(compute_mse(P_GT,P_char));
rmse_ours = sqrt(mean(mse,2));
rmse_best = rmse_ours(pos);

figure(1);
plot(Delta_Set, mean(J,2), 'k-o', 'LineWidth', 2);
xlabel('Threshold \Delta');
ylabel('Cross Validation Cost');
title('NNodes = 200, NBlocks = 10, NObservations = 2, 5 Trials');
grid on;

figure(2);
plot(Delta_Set, rmse_ours, 'k-o', 'LineWidth', 2); hold on;
plot(Delta_Set, rmse_best*ones(length(Delta_Set),1), 'r-x', 'LineWidth', 2); hold on;
plot(Delta_Set, rmse_char*ones(length(Delta_Set),1), 'b-.', 'LineWidth', 2); hold off;
legend('RMSE', 'Cross Validation', 'Chatterjee', 'Location', 'NW');
xlabel('Threshold \Delta');
ylabel('RMSE');
title('NumNodes = 1000, NumBlocks = 10, NumObservations = 2, 50 Trials');
grid on;

