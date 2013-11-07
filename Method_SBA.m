function [H_SBA P_SBA] = Method_SBA(G,Delta)
clusters_SBA     = estimate_blocks_directed(G,Delta);
[H_SBA P_SBA]    = histogram3D(G,clusters_SBA);
