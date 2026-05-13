function L2_norm = compute_vec_L2(e,dx)
% compute_vec_L2: Computes L2 norm of a 3D vector-valued PDE solution over time
% 
% Input:
%   e  - solution matrix of size [3*(N+1) × (M+1)], stacked [u1; u2; u3]
%   dx - spatial discretization step size
%
% Output:
%   L2_norm - vector of L2 norms at each time step

    [rows, Mplus1] = size(e);
    Nplus1 = rows / 3;

    L2_norm = zeros(1, Mplus1);

    for k = 1:Mplus1
        e1 = e(1:3:end, k);
        e2 = e(2:3:end, k);
        e3 = e(3:3:end, k);
        L2_norm(k) = sqrt(sum(e1.^2 + e2.^2 + e3.^2) * dx);
    end
end
