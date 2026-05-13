% Run this file, to have phase portrait for dirichlet and mixed conditions.
% and the error dynamics plot for dirichlet and mixed case. 

% We first define the system parameters. 
a= 1; % Diffusion coefficient
L=1; % Lipschitz constant is always greater than a, then only system will become unstable...
% since eigen values for linearized system are -n^2*pi^2*a/l^2. Here l = pi;
f=@(t,z) L*sin(z);  % nonlinearity function. 
K= 10; % Global control gain.
k=15;  % First Local control gain which apperars because of boundary agents behaving as leaders.
sigma = 1; % Second Local control gain which apperars because of boundary agents behaving as followers.

% gamma_in represents the curve where agents are located initially.
gamma_in=@(x) [0.5*sin(x);0.5*cos(x) ; 3* ones(size(x))];
% gamma represents the curve where agents are required to be deployed.
gamma = @(x) [0.8*(sin(x)).^3; 
              0.05*(13*cos(x)-5*cos(2*x) -2*cos(3*x) -cos(4*x)); 
              5 * ones(size(x))];  % Filling the third row with constant value 5
tf=20;	% simulation time (see ode45 command)
n=size(gamma(0),1); % Agent dimension 

% Evaluating the solution for error dynamics ode system under mixed
% conditions with linear communication between agents.

%[T_M,Y_M]=soln_ODEs_Mixed_linear(a,f,K,k, sigma,gamma,gamma_in,tf,10); % For 11 agents.

% Evaluating the solution for error dynamics ode system under dirichlet
% conditions with linear communication between agents i.e sigma =0.

%[T_D,Y_D]= soln_ODEs_Mixed_linear(a,f,K,k,0,gamma,gamma_in,tf,10); % For 11 agents

% Plotting phase portrait
%figure('name','Phase portrait under dirichlet conditions.')
 %phase_Plot2(Y_D,gamma,gamma_in)
% figure('name','Phase portrait under wentzell conditions.')
 % phase_Plot2(Y_M,gamma,gamma_in)

% Plotting the error graph
figure('name',' Square of the l2-norm of the error')
[T_M_40,Y_M_40]=soln_ODEs_Mixed_linear(a,f,K,k, sigma,gamma,gamma_in,tf,42); % for 40 agents (Linear wentzell)
[T_D_40,Y_D_40]=soln_ODEs_Mixed_linear(a,f,K,k, 0,gamma,gamma_in,tf,42); % for 40 agents (Linear dirichlet )
[CT_M_40,CY_M_40]=soln_ODEs_Mixed_constant(a,f,K,k, sigma,gamma,gamma_in,tf,42); % for 40 agents (Constant wentzell)
[CT_D_40,CY_D_40]=soln_ODEs_Mixed_constant(a,f,K,k, 0,gamma,gamma_in,tf,42); % for 40 agents (Constant Dirichlet)

N=42; % Error plot for 41 agents

% Define colors
black = [0 0 0];
blue = [0 0 1];
green = [0 0.5 0];
red = [1 0 0];

% Plot the data
plot(T_M_40, sum(Y_M_40.^2,2)/(N+1), 'Color', black,'LineWidth', 2);
hold on;
plot(T_D_40, sum(Y_D_40.^2,2)/(N+1), 'Color', green, 'LineWidth', 2);
plot(CT_M_40, sum(CY_M_40.^2,2)/(N+1), 'Color', blue, 'LineWidth', 2);
plot(CT_D_40, sum(CY_D_40.^2,2)/(N+1), 'Color', red,'LineWidth', 2);
hold off;
xlabel('$t$','Interpreter','latex', 'FontSize',25), ylabel('$\|e\|^2$','Interpreter','latex', 'FontSize',25) 
legend('WBC (Communicating)', ...
       'Dirichlet (Communicating)', ...
       'WBC (Isoated)', ...
       'Dirichlet (Isolated)', ...
       'FontSize', 10);
