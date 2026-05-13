% Run this file, to have phase portrait for dirichlet and mixed conditions.
% and the error dynamics plot for dirichlet and mixed case. 

% We first define the system parameters. 
a= 0.003; % Diffusion coefficient
L=0.3; % Lipschitz constant is always greater than a, then only system will become unstable...
% since eigen values for linearized system are -n^2*pi^2*a/l^2. Here l = pi;
f=@(t,z) L*sin(z);  % nonlinearity function. 
K= 1.4; % Global control gain.
k=2;  % First Local control gain which apperars because of boundary agents behaving as leaders.
sigma = 0.8; % Second Local control gain which apperars because of boundary agents behaving as followers.

% gamma_in represents the curve where agents are located initially.
%gamma_in=@(x) [sin(x);cos(x) ; 3* ones(size(x))];
gamma_in=@(x) [0*x;0*x ; 0*x];
% gamma represents the curve where agents are required to be deployed.
gamma = @(x) [1.5*sin(x); 1.5*cos(x); 6*ones(size(x))];

%gamma = @(x) [0.8*(sin(x)).^3; 
             % 0.05*(13*cos(x)-5*cos(2*x) -2*cos(3*x) -cos(4*x)); 
             % 5 * ones(size(x))];  % Filling the third row with constant value 5
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
plot(T_M_40, sum(Y_M_40.^2,2)/(N+1), 'Color', black,'LineWidth', 1);
hold on;
%plot(T_D_40, sum(Y_D_40.^2,2)/(N+1), 'Color', green, 'LineWidth', 3);
plot(CT_M_40, sum(CY_M_40.^2,2)/(N+1), 'Color', blue, 'LineWidth', 1);
%plot(CT_D_40, sum(CY_D_40.^2,2)/(N+1), 'Color', red,'LineWidth', 3);
hold off;




% hold off; % Disable hold on for main plot

xlabel('$t$','Interpreter','latex', 'FontSize',35), ylabel('$\|e(t, \cdot)\|^2$','Interpreter','latex', 'FontSize',35) 
legend('WBC (Communicating)', ...
       'Dirichlet (Communicating)', ...
       'WBC (Isoated)', ...
       'Dirichlet (Isolated)', ...
       'FontSize', 10);
% pHASE pORTRAIT
figure('name','Phase portrait under wentzell conditions.')
phase_Plot2(Y_M_40,gamma,gamma_in)

%% Graph illustrating the comparison between Wentzell boundaries with and without communication across various sigma values.
figure('name',' Variations in error across various sigma values.')
% [LTs0,LEs0]=soln_ODEs_Mixed_linear(a,f,K,k, 0,gamma,gamma_in,tf,42);
[LTs1,LEs1]=soln_ODEs_Mixed_linear(a,f,K,k, 0.1,gamma,gamma_in,tf,42); % for 40 agents (Linear wentzell) sigma = 0.1
[LTs2, LEs2]=soln_ODEs_Mixed_linear(a,f,K,k, 0.4,gamma,gamma_in,tf,42); %  sigma = 0.4
[LTs3, LEs3]=soln_ODEs_Mixed_linear(a,f,K,k, 1.2,gamma,gamma_in,tf,42); %  sigma = 1.2
[LTs4, LEs4]=soln_ODEs_Mixed_linear(a,f,K,k, 100,gamma,gamma_in,tf,42); %  sigma = 100

[CTs1,CEs1]=soln_ODEs_Mixed_constant(a,f,K,k, 0.1,gamma,gamma_in,tf,42); % for 40 agents (Constant wentzell) sigma = 0.1
[CTs2,CEs2]=soln_ODEs_Mixed_constant(a,f,K,k, 0.4,gamma,gamma_in,tf,42); % sigma = 0.4
[CTs3,CEs3]=soln_ODEs_Mixed_constant(a,f,K,k, 1.2,gamma,gamma_in,tf,42); % sigma = 1.2
[CTs4,CEs4]=soln_ODEs_Mixed_constant(a,f,K,k, 100,gamma,gamma_in,tf,42); % sigma = 100
N = 40;
% Plot the data
figure;
plot(LTs1, sum(LEs1.^2, 2)/(N+1), 'Color', 'black', 'LineWidth', 3);
hold on;
plot(LTs2, sum(LEs2.^2, 2)/(N+1), 'Color', 'green', 'LineWidth', 3);
plot(LTs3, sum(LEs3.^2, 2)/(N+1), 'Color', 'blue', 'LineWidth', 3);
plot(LTs4, sum(LEs4.^2, 2)/(N+1), 'Color', 'red', 'LineWidth', 3);
% plot(LTs0, sum(LEs0.^2,2)/(N+1), 'Color', 'orange', 'LineWidth', 3);

plot(CTs1, sum(CEs1.^2, 2)/(N+1), '--', 'Color', 'black', 'LineWidth', 3);
plot(CTs2, sum(CEs2.^2, 2)/(N+1), '--', 'Color', 'green', 'LineWidth', 3);
plot(CTs3, sum(CEs3.^2, 2)/(N+1), '--', 'Color', 'blue', 'LineWidth', 3);
plot(CTs4, sum(CEs4.^2, 2)/(N+1), '--', 'Color', 'red', 'LineWidth', 3);

hold off;

% Setting axes properties
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('$\|e(t, \cdot)\|^2$', 'Interpreter', 'latex', 'FontSize', 30);
set(gca, 'FontSize', 25, 'LineWidth', 3);

% Setting legend properties with two columns
legend({'Interacting, \sigma = 0.1', 'Interacting, \sigma = 0.4', 'Interacting, \sigma = 1.2', 'Interacting, \sigma = 100',...
    'Isolated, \sigma = 0.1', 'Isolated, \sigma = 0.4', 'Isolated, \sigma = 1.2', 'Isolated, \sigma = 100'}, ...
    'Location', 'best', 'FontSize', 28, 'NumColumns', 2);

% Adding box with line width 3
box on;

% Zooming in on the range from t = 1.4 to t = 2
xlim([1.4, 2]);

% Ensuring the figure is displayed correctly
%grid on;







