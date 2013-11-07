function [G P u] = construct_a_graph(w,n,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [G P] construct_a_graph(w,n,T)
% constructs T random graphs G of n nodes based on the stochastic blockmodel w
% 
% Input: w - block model, and KxK matrix
%        n - number of nodes
%        T - number of observations
%
% Output: G - random graph, in dimension n x n x T
%         P - probability of generating an edge; size n x n
%             (needed for computing mse)
% 
% Stanley Chan @ Harvard
% Feb 12, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% % For debugging ...
% clear all
% close all
% clc
% 
% w = [0.9 0.2;
%     0.2 0.7];
% n = 200;
% T = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct a sequence of n nodes
K = size(w,1);                % number of blocks
v = sort(rand(n,1));
u = round((K-1)*v)+1; % sample n ui's from U(0,1)
                              % round off to K-1
                              
% Construct the Probability Matrix (nxn)
P = zeros(n,n);
for i=1:n
    for j=1:n
        ui = u(i);
        uj = u(j);
%         if i==j
%             P(i,j) = 0;
%         else
            P(i,j) = w(ui,uj);
%         end
    end
end


% Construct a random graph with T observations
G = zeros(n,n,T);
for t=1:T
    G(:,:,t) = rand(n,n)<P;
end