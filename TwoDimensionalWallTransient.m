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
dimension = [5 7]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes - 1);

% mesh Fourier number
tau = meshFourier(alpha, dt, dx, dy); 
criteria = 1-4*tau;

%% Explicit Stability Criterion
if  criteria < 0
    warning("Will not converge, consider decreasing dt")
end 


TwoDNodes = NodeSystem(xNodes, yNodes);

%% Referring to Node system

% Upper Left corner = TwoDNodes(1)
% Bottom Left corner = TwoDNodes(end,1)
% Upper Right corner = TwoDNodes(1,end)
% Bottom Right corner = TwoDNodes(end,1)

% Top = TwoDNodes(1,2:yNodes-1)
% Bottom = TwoDNodes(end,2:yNodes-1)
% Left side = TwoDNodes(2:xNodes-1,1)
% Right Side = TwoDNodes(2:xNodes-1,end)

function alpha = ThermalDiffusivity(rho, cp, k)
alpha = k/(rho*cp);
end 

function tau = meshFourier(alpha, dt, dx, dy)
tau = (alpha*dt)/(dx*dy);
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
