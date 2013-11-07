function mse = compute_mse(P,Q)
mse = norm(P(:)-Q(:))^2/numel(P);