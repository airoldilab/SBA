function B = estimate_blocks_largest_gap(G,Q)
%%%%%%%%%%%%%%%%%%%%%%
% Classification and estimation in the SBM based on empirical degrees
% A. Channarond, J. Daudin and S. Robin
% Electronic Journal of Statistics, Vol.6, pp.2574-2601, 2012
%
% Stanley Chan @ Harvard
% May 6, 2012
%
% Input: G - graph
%        Q - number of blocks
%%%%%%%%%%%%%%%%%%%%%%
if size(G,3)>1
    G = sum(G,3);
end
n = size(G,1);

Deg           = sum(G-diag(diag(G)),2);       % degree
DegNorm       = Deg/(n-1);                    % normalized degree
[DegNorm idx] = sort(DegNorm);                % sort normalized degree

DegNormDiff   = diff(DegNorm);                % difference
[~, I]         = sort(DegNormDiff,'descend'); % find Q-1 largest gaps
ii            = [0 sort(I(1:Q-1))' n];               % ii = Q-1 indices of largest gaps

B = cell(Q,1);
for q=1:Q
    B{q} = idx(ii(q)+1:ii(q+1));
end