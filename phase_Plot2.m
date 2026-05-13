function phase_Plot2(y5,gamma, gamma_in)
% Inputs needed are a) soluton of ODE, b) Initial position of agents-curve
% i.e gamma_in c) Target position of agents curve i.e. gamma




N=42; % the number of agents; % the number of agents
% y are the errors, x1 and x2 are the original coordinates 
h = 1/(N);
gammaVal=gamma((0:N)*h);% target curve 
x1=gammaVal(1,:) + y5(:,1:3:end); 
x2=gammaVal(2,:) + y5(:,2:3:end); 
x3=gammaVal(3,:) + y5(:,3:3:end);
gammaPlot=gamma(0:.01:1); 
gammainPlot=gamma_in(0:.01:1);

% Get the size of the matrices
[~, numColumns] = size(x1);


% Define the indices of the columns to be plotted in red
redColumns = [1, 8,15,22,29,36, numColumns];

% Calculate the index for the midpoint
plot3(gammaPlot(1,:),gammaPlot(2,:),gammaPlot(3,:),'Color', [1 0.5 0],'LineWidth',5)  % plot the target curve 
hold on
plot3(gammainPlot(1,:),gammainPlot(2,:),gammainPlot(3,:),'Color', [0 0.5 0],'LineWidth',2)  % plot the initial  curve 
hold on
for i = redColumns
   plot3(x1(:,i), x2(:,i), x3(:,i), 'g', 'LineWidth', 2);
end

% Plot the remaining columns in black
for i = 1:numColumns
    if ~ismember(i, redColumns)
        plot3(x1(:,i), x2(:,i), x3(:,i), 'k');
    end
end
hold off;
% legend('Target Curve', 'Initial Curve', 'Trajectories');
xlabel('$z_i^1$','Interpreter','latex','FontSize',20),...
    ylabel('$z_i^2$','Interpreter','latex','FontSize',20),...
    zlabel('$z_i^3$','Interpreter','latex','FontSize',20)














