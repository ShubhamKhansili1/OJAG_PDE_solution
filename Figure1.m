function Fig1 = Figure1(N, M, T, a1, a2, kappa, K, sigma, num_leaders, gamma, gamma_in,f)
%RUNANDPLOTERROR  Solve the OJA–G parabolic PDEs and plot log-error curves.
%
%   h = RUNANDPLOTERROR(N,M,T,a1,a2,kappa,K,sigma,num_leaders) carries out
%   the following steps for the two diffusion coefficients a1 and a2:
%       1. Calls OJAG_solveParabolicPDE_vec      (with communication)
%       2. Calls OJAG_solveParabolicPDE_constant_vec (no communication)
%       3. Repeats 1 & 2 with σ forced to zero
%       4. Computes the vector-valued L2 error norm  ||e(t,·)||   for each
%          solution, converts it to log-scale, and plots the results in a
%          two-row figure (upper row: a1, lower row: a2).
%
%   INPUTS
%   ──────────────────────────────────────────────────────────────────────
%   N, M, T          : Grid points, time steps, final time (as in your code)
%   a1, a2           : Two diffusion parameters you want to compare
%   kappa, K, sigma  : Control gains and noise parameter
%   num_leaders      : Number of leader agents
%
%   OUTPUT
%   hFig (optional)  : Handle to the generated figure
%
%   Example
%   ──────────────────────────────────────────────────────────────────────
%   N  = 71;   M = 2000;   T = 4;
%   a1 = 0.004; a2 = 0.006;
%   kappa = 2; K = 1.4; sigma = 0.5;   num_leaders = 7;
%
%   runAndPlotError(N,M,T,a1,a2,kappa,K,sigma,num_leaders);
%
%   Author:  <your-name>   •  <date>
% ────────────────────────────────────────────────────────────────────────

%–––– 1.  Preparations –––––––––––––––––––––––––––––––––––––––––––––––––––
dx   = 1/N;

% Anonymous helper that converts the agent-major state matrix E
% (size 3(N+1) × Nt) into ||e(t,·)|| for each time index t:
compute_vec_L2 = @(E) arrayfun(@(tt) ...
    sqrt( ...
        sum( arrayfun(@(i) ...
            norm(E((i-1)*3+1:i*3,tt))^2, 1:N+1 ) ) * dx ), ...
        1:size(E,2) );

%–––– 2.  Solve for a = a1 –––––––––––––––––––––––––––––––––––––––––––––––
[t1 , e_L ]  = OJAG_solveParabolicPDE_vec         (N,M,T,a1,kappa,K,sigma, ...
                                                   gamma,gamma_in,f,num_leaders);
[t2 , e_C ]  = OJAG_solveParabolicPDE_constant_vec(N,M,T,a1,kappa,K,sigma, ...
                                                   gamma,gamma_in,f,num_leaders);
[t11, e_L0]  = OJAG_solveParabolicPDE_vec         (N,M,T,a1,kappa,K,0, ...
                                                   gamma,gamma_in,f,num_leaders);
[t22, e_C0]  = OJAG_solveParabolicPDE_constant_vec(N,M,T,a1,kappa,K,0, ...
                                                   gamma,gamma_in,f,num_leaders);

%–––– 3.  Solve for a = a2 –––––––––––––––––––––––––––––––––––––––––––––––
[t1A , e_LA ] = OJAG_solveParabolicPDE_vec         (N,M,T,a2,kappa,K,sigma, ...
                                                   gamma,gamma_in,f,num_leaders);
[t2A , e_CA ] = OJAG_solveParabolicPDE_constant_vec(N,M,T,a2,kappa,K,sigma, ...
                                                   gamma,gamma_in,f,num_leaders);
[t11A, e_LA0] = OJAG_solveParabolicPDE_vec         (N,M,T,a2,kappa,K,0, ...
                                                   gamma,gamma_in,f,num_leaders);
[t22A, e_CA0] = OJAG_solveParabolicPDE_constant_vec(N,M,T,a2,kappa,K,0, ...
                                                   gamma,gamma_in,f,num_leaders);

%–––– 4.  Post-process (L2 norms & logs) –––––––––––––––––––––––––––––––––
log_L2_L   = log(compute_vec_L2(e_L ));
log_L2_C   = log(compute_vec_L2(e_C ));
log_L2_L0  = log(compute_vec_L2(e_L0));
log_L2_C0  = log(compute_vec_L2(e_C0));

log_L2_LA  = log(compute_vec_L2(e_LA ));
log_L2_CA  = log(compute_vec_L2(e_CA ));
log_L2_LA0 = log(compute_vec_L2(e_LA0));
log_L2_CA0 = log(compute_vec_L2(e_CA0));

%–––– 5.  Plot : All Four Strategy Together–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
Fig1 = figure('Color','w');

% Upper subplot (a = a1)
subplot(2,1,1); hold on; grid on;
plot(t1 ,  log_L2_L  ,'Color',[1 .5 0],'LineWidth',2);
plot(t2 ,  log_L2_C  ,'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot(t11,  log_L2_L0 ,'r','LineWidth',2);
plot(t22,  log_L2_C0 ,'r','LineStyle','--','LineWidth',2);

xlabel('$t$','Interpreter','latex');
ylabel('$\log\!\bigl(\|e(t,\cdot)\|\bigr)$','Interpreter','latex');
title(sprintf('$\\log\\bigl(\\|e(t,\\cdot)\\|\\bigr)$ vs. $t$: $a = %.3f$',a1), ...
      'Interpreter','latex');
legend({'$\sigma \neq 0$ with comm.', '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.',   '$\sigma = 0$ no comm.'}, ...
        'Interpreter','latex','Location','best','FontSize',12);

% Lower subplot (a = a2)
subplot(2,1,2); hold on; grid on;
plot(t1A , log_L2_LA ,'Color',[1 .5 0],'LineWidth',2);
plot(t2A , log_L2_CA ,'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot(t11A, log_L2_LA0,'r','LineWidth',2);
plot(t22A, log_L2_CA0,'r','LineStyle','--','LineWidth',2);

xlabel('$t$','Interpreter','latex');
ylabel('$\log\!\bigl(\|e(t,\cdot)\|\bigr)$','Interpreter','latex');
title(sprintf('$\\log\\bigl(\\|e(t,\\cdot)\\|\\bigr)$ vs. $t$: $a = %.3f$',a2), ...
      'Interpreter','latex');
legend({'$\sigma \neq 0$ with comm.', '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.',   '$\sigma = 0$ no comm.'}, ...
        'Interpreter','latex','Location','best','FontSize',12);

end

%–––– 6.  Plot dirichlet com vs non com, and wentzell com versus non com –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
Fig1 = figure('Color','w');

% Upper subplot (a = a1)
subplot(2,1,1); hold on; grid on;
plot(t1 ,  log_L2_L  ,'Color',[1 .5 0],'LineWidth',2);
plot(t2 ,  log_L2_C  ,'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot(t11,  log_L2_L0 ,'r','LineWidth',2);
plot(t22,  log_L2_C0 ,'r','LineStyle','--','LineWidth',2);

xlabel('$t$','Interpreter','latex');
ylabel('$\log\!\bigl(\|e(t,\cdot)\|\bigr)$','Interpreter','latex');
title(sprintf('$\\log\\bigl(\\|e(t,\\cdot)\\|\\bigr)$ vs. $t$: $a = %.3f$',a1), ...
      'Interpreter','latex');
legend({'$\sigma \neq 0$ with comm.', '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.',   '$\sigma = 0$ no comm.'}, ...
        'Interpreter','latex','Location','best','FontSize',12);

% Lower subplot (a = a2)
subplot(2,1,2); hold on; grid on;
plot(t1A , log_L2_LA ,'Color',[1 .5 0],'LineWidth',2);
plot(t2A , log_L2_CA ,'Color',[1 .5 0],'LineStyle','--','LineWidth',2);
plot(t11A, log_L2_LA0,'r','LineWidth',2);
plot(t22A, log_L2_CA0,'r','LineStyle','--','LineWidth',2);

xlabel('$t$','Interpreter','latex');
ylabel('$\log\!\bigl(\|e(t,\cdot)\|\bigr)$','Interpreter','latex');
title(sprintf('$\\log\\bigl(\\|e(t,\\cdot)\\|\\bigr)$ vs. $t$: $a = %.3f$',a2), ...
      'Interpreter','latex');
legend({'$\sigma \neq 0$ with comm.', '$\sigma \neq 0$ no comm.', ...
        '$\sigma = 0$ with comm.',   '$\sigma = 0$ no comm.'}, ...
        'Interpreter','latex','Location','best','FontSize',12);

end
