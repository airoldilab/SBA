function G = construct_a_graph_from_P(P,n,T)
G = zeros(n,n,T);
for t=1:T
    G(:,:,t) = rand(n,n)<P;
end