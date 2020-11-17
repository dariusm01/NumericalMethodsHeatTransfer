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
dt = 0.0000001; % size of steps
timeSteps = 50000; % number of steps

%% Nodes (horizontal & vertical)
dimension = [3 3]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

xNodes = dimension(2); % Across
dx = L/(xNodes - 1);

yNodes = dimension(1); % Down
dy = H/(yNodes - 1);

rows = yNodes;
cols = xNodes;


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
    temps(1,1) = (4*dt*(0.2500*dx*dy*egen + 0.5000*dx*h3*(- Ts(1,1) + Tinf3) +...
                 0.5000*dy*h1*(- Ts(1,1) + Tinf1) + (0.5000*dx*k*(Ts(2,1) - T(1,1)))/dy -...
                (0.5000*dy*k*(- Ts(1,2) + Ts(1,1)))/dx + (0.2500*Ts(1,1)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

    % Bottom Left corner 
    temps(end,1) = (4*dt*(0.2500*dx*dy*egen + 0.5000*dx*h4*(- Ts(end,1) + Tinf4) +...
                    0.5000*dy*h1*(- Ts(end,1) + Tinf1) - (0.5000*dx*k*(- Ts(end-1,1) +...
                    Ts(end,1)))/dy - (0.5000*dy*k*(- Ts(end,2) + Ts(end,1)))/dx +...
                    (0.2500*Ts(end,1)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

    % Upper Right corner
    temps(1,end) = (4*dt*(0.2500*dx*dy*egen + 0.5000*dx*h3*(- Ts(1,end) +...
                    Tinf3) + 0.5000*dy*h2*(- Ts(1,end) + Tinf2) +...
                    (0.5000*dy*k*(Ts(1,end-1) - Ts(1,end)))/dx +...
                    (0.5000*dx*k*(Ts(2,end) - Ts(1,end)))/dy +...
                    (0.2500*Ts(1,end)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

    % Bottom Right corner 
    temps(end,end) = (4*dt*(0.2500*dx*dy*egen +...
                     0.5000*dx*h4*(- Ts(end,end) + Tinf4) +...
                     0.5000*dy*h2*(- Ts(end,end) + Tinf2) +...
                     (0.5000*dy*k*(Ts(end,end-1) - Ts(end,end)))/dx -...
                     (0.5000*dx*k*(- Ts(end-1,end) +...
                     Ts(end,end)))/dy +...
                     (0.2500*Ts(end,end)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);


    %% Populating the edges
    for i= 2:rows-1
        for j = 2:cols-1

            % Top
            temps(1,j) = (2*dt*(0.5000*dx*dy*egen + 0.5000*dx*h3*(- Ts(1,j) + Tinf3) +...
                         (dy*k*(Ts(1,j-1) - Ts(1,j)))/dx - (dx*k*(- Ts(i+1,j) +...
                         Ts(1,j)))/dy - (dy*k*(- Ts(1,j+1)+ T(1,j)))/dx +...
                         (0.5000*Ts(1,j)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

            % Bottom    
            temps(end,j) = (2*dt*(0.5000*dx*dy*egen + 0.5000*dx*h4*(- Ts(end,j) + Tinf4) +...
                           (dy*k*(Ts(end,j-1) - Ts(end,j)))/dx - (dx*k*(- Ts(end-1,j) +...
                           Ts(end,j)))/dy - (dy*k*(- Ts(end,j+1) + Ts(end,j)))/dx +...
                           (0.5000*Ts(end,j)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

            % Left side
            temps(i,1) = (2*dt*(0.5000*dx*dy*egen + 0.5000*dy*h1*(- Ts(i,1) + Tinf1) - (dx*k*(- Ts(i-1,1) +...
                         Ts(i,1)))/dy - (dx*k*(- Ts(i+1,1) + Ts(i,1)))/dy - (dy*k*(- Ts(i,2) +...
                         Ts(i,1)))/dx + (0.5000*Ts(i,1)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);

            % Right side
            temps(i,end) = (2*dt*(0.5000*dx*dy*egen + 0.5000*dy*h2*(- Ts(i,end) + Tinf2) +...
                           (dy*k*(Ts(i,end-1) - Ts(i,end)))/dx - (dx*k*(- Ts(i-1,end) +...
                           Ts(i,end)))/dy - (dx*k*(- Ts(i+1,end) + Ts(i,end)))/dy +...
                           (0.5000*Ts(i,end)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);   

          %% Interior Nodes
            temps(i,j) = (dt*(dx*dy*egen + (dy*k*(Ts(i,j-1) - Ts(i,j)))/dx +...
                         (dx*k*(Ts(i+1,j) - Ts(i,j)))/dy - (dx*k*(- Ts(i-1,j) +...
                         Ts(i,j)))/dy - (dy*k*(- T(i,j+1) +...
                         Ts(i,j)))/dx + (Ts(i,j)*cp*dx*dy*rho)/dt))/(cp*dx*dy*rho);
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










