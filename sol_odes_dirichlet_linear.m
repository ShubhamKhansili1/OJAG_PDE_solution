function [tmesh,y]=sol_odes_dirichlet_linear(a,f,K,k,gamma,gamma_in,tf,N)
% ODE modelling of MAS with a delay in the leader 
disp('Solving ODEs...')
tic

% See equation (7) in paper, where we defined \(F(t,e)\) =\(f(t, \gamma(x) + e) - f(t, \gamma(x))\),
% Therefore, we have
deltaf=@(t,gamma,y) f(t,gamma+y)-f(t,gamma); 


% Lets solve the dynamics of the first and last agent using ODE45
first_leader=ode45(@(t,y) -k*y+deltaf(t,gamma(0),y),[0,tf],gamma_in(0)-gamma(0)); 

last_leader=ode45(@(t,y) -k*y+deltaf(t,gamma(pi),y),[0,tf],gamma_in(pi)-gamma(pi)); 

% Here we considered the spatial domain as [0, pi]. Therefore we take h as pi/N.
h=pi/(N); 

% we respresent the ode system as 
% y' = Ty + F
% Here, T is laplacian matrix.
T= a/(h^2)*(-2*diag(ones(1,N-1))+diag(ones(1,N-2),-1)+ diag(ones(1,N-2),1));

% lets define the dimension of the state as I_n. Here we took it as 3.
I_n=eye(size(gamma(0),1)); 

% Target position for the agents are defined as 
gammaVal=gamma((1:N-1)*h); 

% Initial position for the agents 
gammainVal=gamma_in((1:N-1)*h); 

% Defining the linear connection between the agents.
% The linear function for x in domain [a,b] is given by
% f(x) = f(a)(b-x/b-a) + f(b)(x-a/b-a). 

% Lets define the linear function for x in [0,pi/2]
    b0= @(x) -(x-0.5*pi)/(0.5*pi);
    b1= @(x) (x)/(0.5*pi);

% Lets define the linear function for x in [pi/2,pi]

     b2 =@(x) -(x-pi)/(0.5*pi);
     b3 =@(x) (x-0.5*pi)/(0.5*pi);
% Creating array which contains these function values
b0val = b0((1: (N/2))*h) ;
b1val = b1((1: (N/2))*h) ;
b2val = b2(((N/2)+1: (N-1))*h) ;
b3val = b3(((N/2)+1: (N-1))*h) ;
M = zeros(N-1);
M(:, N/2) = [b1val';b2val'];

% Using ode45 command.

rhs=@(t,y) kron(T,I_n)*y+deltaf(t,gammaVal(:),y)+kron([a/(h^2);zeros(N-2,1)],I_n)*deval(first_leader,t) + kron([zeros(N-2,1);a/(h^2)],I_n)*deval(last_leader,t)...
    -K*kron([b0val';zeros(N/2-1,1)],I_n)*deval(first_leader,t) -K*kron(M,I_n)*y  -K*kron([zeros((N/2),1);b3val'],I_n)*deval(last_leader,t); 
[tmesh,y]=ode45(rhs,[0,tf], gammainVal(:)-gammaVal(:)); 
%% Combining data 
y=[deval(first_leader,tmesh)',y, deval(last_leader,tmesh)']; 
disp("ODEs are solved in " + num2str(toc) + " sec")