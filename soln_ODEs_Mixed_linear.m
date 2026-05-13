function [tmesh,y]=soln_ODEs_Mixed_linear(a,f,K,k, sigma,gamma,gamma_in,tf,N)
% This functon file solves the ODE system (error dynamics (denoted as y here), see paper) using
% ode45 command. 

disp('Solving ODEs...')
tic
% See equation (7) in paper, where we defined \(F(t,e)\) =\(f(t, \gamma(x) + e) - f(t, \gamma(x))\),
% Therefore, we have
deltaf=@(t,gamma,y) f(t,gamma+y)-f(t,gamma); 

% Here we considered the spatial domain as [0, pi]. Therefore we take h as pi/N.
% Change h, in "Laplacian_Matrix_for_mixed_BC" and in "Linear_connection_mixed(N)" if you change domain.

h=1/(N); 

% we respresent the ode system as 
% y' = Ty + F
% Here, T is laplacian matrix, see function file.

T= Laplacian_Matrix_for_mixed_BC(N,k,sigma,a); % Laplacian matrix for N+1 agents.

% lets define the dimension of the state as I_n. Here we took it as 3.

I_n=eye(size(gamma(0),1)); 

% Target position for the agents are defined as 
gammaVal=gamma((0:N)*h); 

% Initial position for the followers 
gammainVal=gamma_in((0:N)*h); 

% Defining the linear connection between the agents. See function file.
M = Linear_connection_mixed(N);

% Using ode45 command.

rhs=@(t,y) kron(T,I_n)*y+ deltaf(t,gammaVal(:),y) -K*kron(M,I_n)*y;
[tmesh,y]=ode45(rhs,[0,tf], gammainVal(:) - gammaVal(:)); 

disp("ODEs are solved in " + num2str(toc) + " sec")