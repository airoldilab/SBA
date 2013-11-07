function B = estimate_blocks_directed(G,Delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function estimate_blocks(G,Delta)
% returns the clusters B
% 
% Input: G, an n x n x T random graph
%        Delta, threshold
% Output: B, the clusters
%
% Stanley Chan @ Harvard
% Apr 23, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Goal: estimate w
n = size(G,1);
T = size(G,3);

% Paremeters
B     = cell(1,1);       % Blocks

% Pick a pivot
PivotIdx = randi(n,1);   % Define pivot #1 arbitrarily
B{1}     = PivotIdx(1);  % Block #1 should contain pivot #1

% Initialize the set of un-assigned nodes
NotAssignedVector           = true(n,1);
NotAssignedVector(PivotIdx) = false;
NotAssignedIdx              = find(NotAssignedVector);

% Loop until 
% (1) All indices have been assigned; 
% (2) All nodes have been scanned
t = 1;
while (~isempty(NotAssignedIdx))&&(t<n)
    
    % Pick arbitrarily an i from un-assigned nodes
    if length(NotAssignedIdx)>1
        i = randsample(NotAssignedIdx,1); % randomly pick one i
    else
        i = NotAssignedIdx;               % last case
    end
    NotAssignedVector(i) = false;         % update NotAssignedVector
    NotAssignedIdx = find(NotAssignedVector);

    dhat = zeros(length(PivotIdx),1);
    % Loop through the pivots
    for j=1:length(PivotIdx)
        % Define the jth pivot
        bj = PivotIdx(j);
        
        % Define the set S (neighborhood for computing dhat)
        SVector         = true(n,1);
        SVector([i,bj]) = false;
        SIdx            = find(SVector);

        % Compute dhat
        Term1 = sum( ((1/floor((T+1)/2))*sum(G(i, SIdx,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(i, SIdx,floor((T+1)/2)+1:T),3))  );
        Term2 = sum( ((1/floor((T+1)/2))*sum(G(bj,SIdx,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(bj,SIdx,floor((T+1)/2)+1:T),3))  );
        Term3 = sum( ((1/floor((T+1)/2))*sum(G(i, SIdx,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(bj,SIdx,floor((T+1)/2)+1:T),3))  );
        Term4 = sum( ((1/floor((T+1)/2))*sum(G(bj,SIdx,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(i, SIdx,floor((T+1)/2)+1:T),3))  );

        Term5 = sum( ((1/floor((T+1)/2))*sum(G(SIdx,i ,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(SIdx,i, floor((T+1)/2)+1:T),3))  );
        Term6 = sum( ((1/floor((T+1)/2))*sum(G(SIdx,bj,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(SIdx,bj,floor((T+1)/2)+1:T),3))  );
        Term7 = sum( ((1/floor((T+1)/2))*sum(G(SIdx,i, 1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(SIdx,bj,floor((T+1)/2)+1:T),3))  );
        Term8 = sum( ((1/floor((T+1)/2))*sum(G(SIdx,bj,1:floor((T+1)/2)),3)).*((1/(T-floor((T+1)/2)))*sum(G(SIdx,i, floor((T+1)/2)+1:T),3))  );
        
        dhatTmp = 0.5*(abs(Term1+Term2-Term3-Term4) + abs(Term5+Term6-Term7-Term8));
        dhat(j) = sqrt(abs(dhatTmp/numel(SIdx)));
    end
    
    % Assign Clusters
    % Look for minimum distance
    [Val Idx] = min(dhat);
    if Val<Delta
        % If min distance < Delta, assign to one of the existing blocks
        B{Idx} = [B{Idx} i];
    else
        % If min distance > Delta, make a new block; Put i as pivot
        B{length(PivotIdx)+1} = i;
        PivotIdx = [PivotIdx i];
    end

    t = t+1;
end