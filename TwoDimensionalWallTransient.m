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
timeSteps = 1000; % number of steps


%% Nodes (horizontal & vertical)
dimension = [10 10]; % any # of nodes (x-direction) & nodes (y-direction)
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


TwoDNodes = NodeSystem(rows, cols);
temps = zeros(size(TwoDNodes));
Ts = Tinitial*ones(size(TwoDNodes)); % Initial


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
        temps(i,1) = ((((alpha*dt)/(dx^2))/dy)*((h1*dx^2*(Tinf1-Ts(i,1))/k)) + (dy*(Ts(i,2)-Ts(i,1)))...
            + dx^2*(Ts(i-1,1)+Ts(i+1,1)-2*Ts(i,1))/2*dy + egen*dx^2*dy/2*k) + Ts(i,1); % generic for now, insert eq     

        % Right side
        temps(i,end) = ((((alpha*dt)/(dx^2))/dy)*((h2*dx^2*(Tinf2-Ts(i,end))/k)) + (dy*(Ts(i,end-1)-Ts(i,end)))...
            + dx^2*(Ts(i-1,end)+Ts(i+1,end)-2*Ts(i,end))/2*dy + egen*dx^2*dy/2*k) + Ts(i,end); % generic for now, insert eq    
       
        %Top
        temps(1,j) = ((((alpha*dt)/(dx^2))/dy)*((h3*dx*dy*(Tinf2-Ts(1,j))/k)) + (dx^2*(Ts(2,j)-Ts(1,j))/dy)...
             + dy*(Ts(1,j-1)+Ts(1,j+1)-2*Ts(1,j))/2 + egen*dx^2*dy/2*k) + Ts(1,j); % generic for now, insert eq 
         
         %Bottom
         temps(end,j) = ((((alpha*dt)/(dx^2))/dy)*((h4*dx*dy*(Tinf2-Ts(end,j))/k)) + (dx^2*(Ts(end-1,j)-Ts(end,j))/dy)...
             + dy*(Ts(end,j-1)+Ts(end,j+1)-2*Ts(end,j))/2 + egen*dx^2*dy/2*k) + Ts(end,j); % generic for now, insert eq
  

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



function alpha = ThermalDiffusivity(rho, cp, k)
alpha = k/(rho*cp);
end 

function x = Interpolation(y2, y1, x2, x1, YourVal)

m = (y2-y1)/(x2-x1);

x = ((YourVal - y1)/m) + x1;
end 

function T = KelvintoC(x)
T = x-273.15;
end 

function y = cm_to_m(x)
y = x/100;
end 

function z = NodeSystem(rows, cols)

    Matrix = zeros(rows,cols);
    
    % Top 
    for i = 1:cols
        Matrix(1,i) = i; 
    end 

    % Bottom 
    for i = 1:cols
        Matrix(end,i) = cols + i;
    end 

    % Left Side
    for i = 2:rows-1
        Matrix(i,1) = 2*cols + (i-1);
    end 

    % Right Side
    for i = 2:rows-1
        Matrix(i,end) = 2*cols + (rows-3) + i;
    end 

    A = 1;
    
    % Inside
    for i=2:rows-1
        for j=2:cols -1
            Matrix(i,j) = (2*cols)+ (2*(rows-2))  + A;
            A = A+1;
        end 
    end 
    
    z = Matrix;
    
end 


function tau = meshFourier(alpha, dt, dx)
tau = (alpha*dt)/(dx^2);
end 







