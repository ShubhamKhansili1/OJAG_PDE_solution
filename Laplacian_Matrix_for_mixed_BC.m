function matrix = Laplacian_Matrix_for_mixed_BC(N,k,sigma,a)

% The function file creates laplacian matrix of size N+1 times N+1 for
% i = 0,1,2,...N. 

% Let the local control gain be k and sigma and h = 1/N (see boundary conditions in paper)
    h=1/N;
% Initialize the matrix with zeros.
    matrix = zeros(N+1);

% Since, the first and last agent satisfy the following BC

% e_{t}(t,j) = -\kappa e(t,j) -(-1)^(j+1) \sigma e_x(t,j) + F(t,e(t,j)),
% ~~t>0, for j =0,1. 

% Therefore, we fill the first and last row of the matrix using above BC.
    matrix(1,1:2) = [-k - (sigma/h), (sigma/h)];
    matrix(N+1,N:N+1) = [(sigma/h), -k - (sigma/h)];
   
% As all intermediate agent follows the equation, a*\Big( \frac{z_{i+1}(t) - 2z_{i}(t)} + z_{i-1}(t)}{h^{2}}\Big)
% Therefore, we fill the middle rows as 
    for i = 2:N
        matrix(i,i-1:i+1) = (a/h^2)*[1, -2, 1];
    end

    % Fill the Nth row
  

    % Display the matrix
end