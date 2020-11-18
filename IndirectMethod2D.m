%% Input paramters
k = 25; % W/mk
egen = 0; % W/m^3
h = 0; % W/m^2k
Tinf = 0; % k
L = 60; % cm 
H = 30; % cm
iter = 0;
dimension = [7 4]; % any # of nodes (x-direction) & nodes (y-direction)

%% How many nodes across and down
xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes -1);

%% Initial Temperature Guesses
T = zeros(xNodes,yNodes).';

%% 2D Boundary Conditions
 
x = 0:dx:xNodes*dx;
x(end) = []; % getting rid of last number

T(1,1:xNodes) = 100*sin((pi*x)/60);  % Top
T(1,end) = 0;
T(end,1:xNodes) = 0; % Bottom
T(1:yNodes,1) = 0; % Left
T(1:yNodes,end) = 0; % right


% Don't do this for the project lmao 
temps = zeros(xNodes,yNodes).';

%% Iterative process
while(1)
    
    T(1,1:xNodes) = 100*sin((pi*x)/60);  % Top
    T(1,end) = 0;
    
    % Interior Nodes (since its a transpose this will be different for the
    % project)
    for m = 2:yNodes-1 % rows
        for n = 2:xNodes-1  % cols
            temps(m,n) = ((k*dy^2*(T(m+1,n)+T(m-1,n))) + (k*dx^2*(T(m,n+1)+T(m,n-1))) +...
                (egen*(dx^2)*(dy^2)))/(k*(2*dy^2 +2*dx^2));
        end 
    end
    
    PercentChange  = (abs((temps(m,n) - T(m,n))/T(m,n)))*100;
    
    if PercentChange <= 0.5e-7 % arbitrary value (just small enough for convergence)  
        break
    end 
    
    % redefining for next iteration
    T = temps;
    iter = iter + 1; 
end

% Won't need this for the project 
temps = flipud(T); % Flipping the top and bottom 

%% Plotting
figure(1)
contourf(temps);
colorbar

figure(2)
mesh(temps);