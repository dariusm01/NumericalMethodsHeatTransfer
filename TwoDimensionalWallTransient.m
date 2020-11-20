%% Adding the functions to the filepath
addpath('ProjectFunctions')

%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 2350; % W/m^2k (left)
h2 = 360; % W/m^2k (right)
h3 = 5; % W/m^2k (top) 
h4 = 45; % W/m^2k (bottom)
Tinf1 = KelvintoC(3350); % 3350k to °C
Tinf2 = KelvintoC(295); % 295k to °C
Tinf3 = Tinf2;
Tinf4 = Tinf2;
L = cm_to_m(1); % 1cm to m 
H = cm_to_m(30); % 30cm to m 
rho = 7900; % kg/m^3
cp = Interpolation(300, 200, 477, 402, 295); % J/kg*k
alpha = ThermalDiffusivity(rho, cp, k); % m^2/s
dt = 1e-3; % size of steps
timeSteps = 750000; % number of steps


%% Nodes (horizontal & vertical)
dimension = [3 3]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

xNodes = dimension(2); % Across
dx = L/(xNodes - 1);

yNodes = dimension(1); % Down
dy = H/(yNodes - 1);

rows = yNodes;
cols = xNodes;


%% Explicit Stability Criterion
tauX = meshFourier(alpha, dt, dx);
tauY = meshFourier(alpha, dt, dy);
criteriaX = 1-4*tauX;
criteriaY = 1-4*tauY;

if  criteriaX < 0 || criteriaY < 0
    warning("Will not converge, consider decreasing dt")
end 

%% Using the explicit approach
nodes = rows*cols; % total number of nodes
T = zeros(timeSteps, nodes);
Tinitial =  KelvintoC(295); % 295k to °C
T(1,:) = Tinitial; % setting the first row to the intial temp. These will get updated down the column

%% If you want to see how the nodes are laid out
TwoDNodes = NodeSystem(rows, cols);

%% Creating array for temperatures
temps = zeros(size(TwoDNodes)); % 'current'
Ts = Tinitial*ones(size(TwoDNodes)); % 'old'


for k = 2:timeSteps
    %% Populating the corners
    
    % Upper Left corner 
    temps(1,1) = (2*tauX*(Ts(1,2)-Ts(1,1))) + (2*tauY*(Ts(2,1)-Ts(1,1)))+...
        (((2*h3*dt)/(rho*cp*dy))*(Tinf3-Ts(1,1))) + (((2*h1*dt)/(rho*cp*dx))*(Tinf1-Ts(1,1)))+...
        ((egen*dt)/(rho*cp)) + Ts(1,1); 

    % Bottom Left corner 
    temps(end,1) = (2*tauX*(Ts(end,2)-Ts(end,1))) + (2*tauY*(Ts(end-1,1)-Ts(end,1)))+...
        (((2*h4*dt)/(rho*cp*dy))*(Tinf4-Ts(end,1))) + (((2*h1*dt)/(rho*cp*dx))*(Tinf1-Ts(end,1)))+...
        ((egen*dt)/(rho*cp)) + Ts(end,1);

    % Upper Right corner
    temps(1,end) = (2*tauX*(Ts(1,end-1)-Ts(1,end))) + (2*tauY*(Ts(2,end)-Ts(1,end)))+...
        (((2*h3*dt)/(rho*cp*dy))*(Tinf3-Ts(1,end))) + (((2*h2*dt)/(rho*cp*dx))*(Tinf2-Ts(1,end)))+...
        ((egen*dt)/(rho*cp)) + Ts(1,end);

    % Bottom Right corner 
    temps(end,end) = (2*tauX*(Ts(end,end-1)-Ts(end,end))) + (2*tauY*(Ts(end-1,end)-Ts(end,end)))+...
        (((2*h4*dt)/(rho*cp*dy))*(Tinf4-Ts(end,end))) + (((2*h2*dt)/(rho*cp*dx))*(Tinf2-Ts(end,end)))+...
        ((egen*dt)/(rho*cp)) + Ts(end,end);


    %% Populating the sides
    for i = 2:rows-1
        for j = 2:cols-1

        % Left side
        temps(i,1) =  (2*tauX*(Ts(i,2)-Ts(i,1))) + (tauY*(Ts(i-1,1)-Ts(i,1)))+...
            (tauY*(Ts(i+1,1)-Ts(i,end))) + (((2*h1*dt)/(rho*cp*dx))*(Tinf1-Ts(i,1))) +...
            ((egen*dt)/(rho*cp)) + Ts(i,1);

        % Right side
        temps(i,end) = (2*tauX*(Ts(i,end-1)-Ts(i,end))) + (tauY*(Ts(i-1,end)-Ts(i,end)))+...
            (tauY*(Ts(i+1,end)-Ts(i,end))) + (((2*h2*dt)/(rho*cp*dx))*(Tinf2-Ts(i,end))) +...
            ((egen*dt)/(rho*cp)) + Ts(i,end);
       
        %Top
        temps(1,j) = (2*tauY*(Ts(2,j)-Ts(1,j))) + (tauX*(Ts(1,j-1)-Ts(1,j)))+...
            (tauX*(Ts(1,j+1)-Ts(1,j))) + (((2*h3*dt)/(rho*cp*dx))*(Tinf3-Ts(1,j))) +...
            ((egen*dt)/(rho*cp)) + Ts(1,j);
        
        %Bottom
        temps(end,j) = (2*tauY*(Ts(end-1,j)-Ts(end,j))) + (tauX*(Ts(end,j-1)-Ts(end,j)))+...
            (tauX*(Ts(end,j+1)-Ts(end,j))) + (((2*h4*dt)/(rho*cp*dx))*(Tinf4-Ts(end,j))) +...
            ((egen*dt)/(rho*cp)) + Ts(end,j);

          %% Interior Nodes
            temps(i,j) = (-Ts(i,j)*(2*tauX + 2*tauY - 1)) + tauX*(Ts(i,j-1) + Ts(i,j+1))+...
                tauY*(Ts(i+1,j) + Ts(i-1,j)) + ((alpha*dt*egen)/k);
        end 
    end
    
    Ts = temps;

    B = reshape(temps,[1,nodes]);


    T(k,:) = B;
end

FinalTempsTransient2D = T;

%% Plotting Temperature History of boundaries
figure(1)
% All rows, first column
plot(FinalTempsTransient2D(:,1))
grid on
xlabel("Time steps")
ylabel("Temperature °C")
title("Temperature History of Combustion Chamber Boundary")

figure(2)
% All rows, last column
plot(FinalTempsTransient2D(:,nodes))
grid on
xlabel("Time steps")
ylabel("Temperature °C")
title("Temperature History of Outside Boundary")






