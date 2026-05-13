function matrix = Constant_connection_mixed(N)

% This function file create a matrix which establishes a linear connection
% between the agents.

% Let i = 0, 1, 2,..., N be the agents. (Keep N = 40)

% Here, we consider total of three leaders. Two are at the boundary and one in middle of domain.
% That is, if we assume N to be even, then the third leader agent 
% will be i = (N/2).



% creating the matrix
 matrix = zeros(N+1);
 matrix(:,1) = [ones(4,1);zeros(39,1)]; % containing first leader values
 matrix(:, 8) = [zeros(4,1); ones(7,1); zeros(32,1)];
 matrix(:, 15) = [zeros(11,1);ones(7,1); zeros(25,1)];
 matrix(:, 22) = [zeros(18,1);ones(7,1); zeros(18,1)];% containing mid leader values
 matrix(:, 29) = [zeros(25,1);ones(7,1); zeros(11,1)];
 matrix(:, 36) = [zeros(32,1);ones(7,1); zeros(4,1)];
matrix(:, 43) = [zeros(39,1);ones(4,1)]; % containing last leader values
