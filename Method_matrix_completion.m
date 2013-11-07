function out = Method_matrix_completion(G)
G   = mean(G,3);
tol = 1e-3;
[X, S, Y, ~] = OptSpace(2*(G-0.5),[],20,tol);
out = (X*S*Y'+1)/2;
out(out>1) = 1;
out(out<0) = 0;
