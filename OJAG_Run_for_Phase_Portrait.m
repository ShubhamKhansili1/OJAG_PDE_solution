
% Run this file for phase portrait.
% system parameters are as follows: 
a= 0.01; % Diffusion coefficient
L=2; % Lipschitz constant is always greater than a, then only system will become unstable...
% since eigen values for linearized system are -n^2*pi^2*a/l^2. Here l = pi;
f=@(t,z) L*z;  % nonlinearity function. 
K= 5; % Global control gain.
k=1;  % First Local control gain which apperars because of boundary agents behaving as leaders.
sigma = 6; % Second Local control gain which apperars because of boundary agents behaving as followers.
num_leaders=7;
N=40; % Number of agents are N+1 
M = 1000;
T = 4;
kappa = 1;

% gamma_in: Initial Curve
 %gamma_in=@(x) [sin(x);cos(x) ; 3* ones(size(x))];
 gamma_in=@(x) [0*x;0*x ; 0*x];

% gamma: Target Curve
gamma = @(x) [5*sin(x); 5*cos(x); 5*ones(size(x))];
%gamma = @(x) [1.5*sin(x); 0*x; 0*x];

%gamma = @(x) [0.8*(sin(x)).^3; 
             % 0.05*(13*cos(x)-5*cos(2*x) -2*cos(3*x) -cos(4*x)); 
             % 5 * ones(size(x))];  % Filling the third row with constant value 5

 tf=20;	% Simulation time (see ode45 MATLAB command)
 
 n=size(gamma(0),1); % Agent dimension 

% Evaluating the solution for the ODE-based MAS under with Wentzell
% boundary control and inter-leader communication.

 %[T_M_40,Y_M_40]=OJAG_soln_ODEs_Wentzell(a,f,K,k, sigma,gamma,gamma_in,tf,N, num_leaders); % for 40 agents (Linear wentzell)

[t, e] =OJAG_solveParabolicPDE_vec(N, M, T, a, kappa, K, sigma, gamma, f, num_leaders);

% Phase portrait
figure('name','Phase portrait under wentzell conditions.')
%OJAG_Phase_Plot(Y_M_40,gamma,gamma_in,N)
OJAG_Phase_Plot_vector(e, gamma, gamma_in, N)



