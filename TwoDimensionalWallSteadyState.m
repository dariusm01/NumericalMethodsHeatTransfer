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
iter = 0; % iteration counter
iterLimit = 10000000;

%% Nodes (horizontal & vertical)
dimension = [3 3]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes - 1);

% Transposed so now x = cols & y = rows
% x is across & y is downward
T = zeros(xNodes,yNodes).'; % 'old' temps

temps = zeros(xNodes,yNodes).'; % 'new' temps


while iter < iterLimit
    
    %% Corners 
    
    % Upper Left Corner                                          
    temps(1,1) = ((h1*dx*(dy^2)*Tinf1) + (k*(dy^2)*T(1,2)) +...
        (k*(dx^2)*T(2,1)) + (h3*(dx^2)*dy*Tinf3)+...
        ((egen/2)*dx^2*dy^2))/((h1*dx*(dy^2)) + k*dy^2 + k*dx^2 + (h3*(dx^2)*dy));

    % Upper Right Corner                          
    temps(1,end) = ((h3*dy*(dx^2)*Tinf3) + (k*(dy^2)*T(1,end-1)) +...
        (k*(dx^2)*T(2,end)) + (h2*(dy^2)*dx*Tinf2)+...
        ((egen/2)*dx^2*dy^2))/((h3*dy*(dx^2)) + k*dx^2 + k*dy^2 + (h2*(dy^2)*dx));

    % Lower Left Corner
    temps(end,1) = ((h1*dx*(dy^2)*Tinf1) + (k*(dy^2)*T(end,2)) +...
        (k*(dx^2)*T(end-1,1)) + (h4*(dx^2)*dy*Tinf4)+...
        ((egen/2)*dx^2*dy^2))/((h1*dx*(dy^2)) + k*dx^2 + k*dy^2 + (h4*(dx^2)*dy));

    % Lower Right Corner
    temps(end,end) = ((h4*dy*(dx^2)*Tinf4) + (k*(dy^2)*T(end,end-1)) +...
        (k*(dx^2)*T(end-1,end)) + (h2*(dy^2)*dx*Tinf2)+...
        ((egen/2)*dx^2*dy^2))/((h4*dy*(dx^2)) + k*dx^2 + k*dy^2 + (h2*(dy^2)*dx));

    %% Sides
    for m = 2:yNodes-1 % y is the rows 
        for n = 2:xNodes-1  % x is the cols
            
        % Top
        temps(1,n) = ((h3*(dx^2)*dy*Tinf3) + (2*k*(dy^2)*(T(m,n-1) + T(m,n+1))) +...
            (2*k*(dx^2)*T(m+1,n)) + (egen*(dx^2)*(dy^2)))/((h3*(dx^2)*dy) + (4*k*(dy^2)) +...
            (2*k*(dx^2)));

        % Bottom
        temps(end,n) = ((h4*(dx^2)*dy*Tinf4) + (2*k*(dy^2)*(T(m,n-1) + T(m,n+1))) +...
            (2*k*(dx^2)*T(m-1,n)) + (egen*(dx^2)*(dy^2)))/((h4*(dx^2)*dy) + (4*k*(dy^2)) +...
            (2*k*(dx^2)));

        % Left
        temps(m,1) = ((h1*(dy^2)*dx*Tinf1) + (2*k*(dx^2)*(T(m-1,n) + T(m+1,n))) +...
            (2*k*(dy^2)*T(m,n+1)) + (egen*(dx^2)*(dy^2)))/((h1*(dy^2)*dx) + (4*k*(dx^2)) +...
            (2*k*(dy^2)));

        % Right
        temps(m,end) = ((h2*(dy^2)*dx*Tinf2) + (2*k*(dx^2)*(T(m-1,n) + T(m+1,n))) +...
            (2*k*(dy^2)*T(m,n-1)) + (egen*(dx^2)*(dy^2)))/((h2*(dy^2)*dx) + (4*k*(dx^2)) +...
            (2*k*(dy^2)));
        end 
    end 
    
    %% Interior nodes    
    for i = 2:yNodes-1
        for j = 2:xNodes-1 
        temps(i,j) = ((k*dy^2*(T(i,j+1)+T(i,j-1))) + (k*dx^2*(T(i+1,j)+T(i-1,j))) +...
            (egen*(dx^2)*(dy^2)))/(k*(2*dy^2 +2*dx^2));
        end 
    end 
    
    T = temps;
    iter = iter + 1; 
    
end 

FinalTempsSteadyState2D = temps;

x = 0:dx:L;
y = 0:dy:H;

%% Plotting
figure(1)
contourf(x,y,FinalTempsSteadyState2D);
title("2D Dimensional Steady State Temperatures Contour")
xlabel("Length (m)")
ylabel("Height (m)")
colorbar

figure(2)
mesh(x,y,FinalTempsSteadyState2D);
title("2D Dimensional Steady State Temperatures Mesh")
xlabel("Length (m)")
ylabel("Height (m)")
colorbar

function T = KelvintoC(x)
T = x-273.15;
end 

function x = Interpolation(y2, y1, x2, x1, YourVal)

m = (y2-y1)/(x2-x1);

x = ((YourVal - y1)/m) + x1;
end 

function y = cm_to_m(x)
y = x/100;
end 