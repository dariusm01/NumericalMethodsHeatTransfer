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
dt = 0.005; % size of steps
timeSteps = 50000; % number of steps

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
% if  criteria < 0
%     warning("Will not converge, consider decreasing dt")
% end 


%% Using the explicit approach
n = rows*cols; % total number of nodes
T = zeros(timeSteps, n);
Tinitial = 295; % k
T(1,:) = Tinitial; % setting the first row to the intial temp. These will get updated down the column


TwoDNodes = NodeSystem(rows, cols);
temps = zeros(size(TwoDNodes));
Ts = 0.5*ones(size(TwoDNodes)); % very generic

%% Referring to Node system
% Top = TwoDNodes(1,2:yNodes-1)
% Bottom = TwoDNodes(end,2:yNodes-1)
% Left side = TwoDNodes(2:xNodes-1,1)
% Right Side = TwoDNodes(2:xNodes-1,end)

%% Populating the corners
% Upper Left corner 
temps(1,1) = (alpha*dt/dx^2*dy^2)*(h1*dx^2*((Tinf1-T(1,1))/dy) + k*(T(2,1)-T(1,1)) + h3*dx*(Tinf2-T(1,1))...
    + (k*dx^2*(T(1,2)-T(1,1))/dy^2) + egen*dx^2/2) + T(1,1); % generic for now, insert eq

% Bottom Left corner 
temps(end,1) = (alpha*dt/dx^2*dy^2)*(h1*dx^2*((Tinf1-T(end,1))/dy) + k*(T(end,2)-T(end,1)) + h4*dx*(Tinf2-T(end,1))...
    + (k*dx^2*(T(end-1,1)-T(end,1))/dy^2) + egen*dx^2/2) + T(end,1); % generic for now, insert eq    

% Upper Right corner
temps(1,end) = (alpha*dt/dx^2*dy^2)*(h2*dx^2*((Tinf2-T(1,end))/dy) + k*(T(1,end-1)-T(1,end)) + h3*dx*(Tinf2-T(1,end))...
    + (k*dx^2*(T(2,end)-T(1,end))/dy^2) + egen*dx^2/2) + T(1,end); % generic for now, insert eq    

% Bottom Right corner 
temps(end,end) = (alpha*dt/dx^2*dy^2)*(h2*dx^2*((Tinf2-T(end,end))/dy) + k*(T(end,end-1)-T(end,end)) + h4*dx*(Tinf2-T(1end,end))...
    + (k*dx^2*(T(end-1,end)-T(end,end))/dy^2) + egen*dx^2/2) + T(end,end); % generic for now, insert eq     


%% Populating the top and bottom
for i=1:rows
    for j = 1:cols
        
        % Top
        if TwoDNodes(i,j) > TwoDNodes(1) && TwoDNodes(i,j) < TwoDNodes(1,end)
            temps(i,j) = 5; % generic for now, insert eq     
            
        % Bottom    
        elseif TwoDNodes(i,j) > TwoDNodes(end,1) && TwoDNodes(i,j) < TwoDNodes(end,end)
            temps(i,j) = 7; % generic for now, insert eq           
        end 
    end 
end 

%% Populating the sides
for i = 2:rows-1
    
    % Left side
    temps(i,1) = 9; % generic for now, insert eq     
    
    % Right side
    temps(i,end) = 11; % generic for now, insert eq     
end 

%% Interior Nodes
for i = 2:rows-1
    for j = 2:cols-1
        
        temps(i,j) = (((alpha*dt)/((dx^2)*(dy^2))) * ((-Ts(i,j)*(2*dy^2 + 2*dx^2)) + (dy^2*(Ts(i,j-1)+Ts(i,j+1)))+... 
            (dx^2*(Ts(i-1,j)+Ts(i+1,j))) + (egen*(dx^2 * dy^2)/k))) + Ts(i,j);
    end 
end 

% % I need a way to place the temps from this matrix into the (steps, nodes) matrix so you 
% % can see how each node updates after each timestep 
% for i = 2:timeSteps
%     T(i,:) = 
% end 

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
