%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 360; % W/m^2k (top)
h2 = 2350; % W/m^2k (bottom)
h3 = 8800; % W/m^2k (coolant) 
Tinf1 = KelvintoC(295); % 295k to °C
Tinf2 = KelvintoC(3350); % 3350k to °C
Tinf3 = KelvintoC(183); % 183k to °C
L = cm_to_m(.15); % .15 cm to m 
H = cm_to_m(.35); % .35 cm to m 
iter = 0; % iteration counter
iterLimit = 10000000;

%% Nodes (horizontal & vertical)
dimension = [4 8]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)
% set to 4 x 8 for now to fit my diagram and the given dimensions

xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes - 1);

% Transposed so now x = cols & y = rows
% x is across & y is downward
T = zeros(xNodes,yNodes).'; % 'old' temps

temps = zeros(xNodes,yNodes).'; % 'new' temps

% Format is T(y,x)

while iter < iterLimit
    
    %% Outside Corners
    
    % Upper Left Corner, Node 1                                          
    temps(1,1) = (k*T(1,2)*(dy^2)+(dx^2)*(h1*Tinf1*dy+k*T(2,1))) / ...
        ((dx^2)*(h1*dy+k)+k*(dy^2));

    % Upper Right Corner, Node 2                          
    temps(1,end) = (k*T(1,end-1)*(dy^2)+(dx^2)*(h1*Tinf1*dy+k*T(2,end))) / ...
        ((dx^2)*(h1*dy+k)+k*(dy^2));

    % Lower Left Corner, Node 8
    temps(end,1) = (k*T(end,2)*(dy^2)+(dx^2)*(h2*Tinf2*dy+k*T(end-1,1))) / ...
        ((dx^2)*(h2*dy+k)+k*(dy^2));

    % Lower Right Corner, Node 7
    temps(end,end) = (k*T(end,end-1)*(dy^2)+(dx^2)*(h2*Tinf2*dy+k*T(end-1,end))) / ...
        ((dx^2)*(h2*dy+k)+k*(dy^2));

    %% Inside Corners
    % T(y,x) = T(m,n), so m = y and n = x
    % Upper Left Corner, Node 4
    temps(3,2) = (3*k*T(2,2)*(dy^2)+dx*(dx*(h3*Tinf3*dy+k*(2*T(2,2)*T(4,2))+(h3*Tinf3*(dy^2))))) / ...
        ((dx^2)*(h3*dy+3*k)+h3*dx*(dy^2)+3*k*(dy^2));
    
    % Upper Right Corner, Node 3
    temps(3,4) = (k*T(3,3)*(dy^2)+(dx^2)*(h3*Tinf3*dy+k*T(2,4))) / ...
        ((dx^2)*(h3*dy+k)+k*(dy^2));
    
    % Lower Left Corner, Node 5
    temps(7,2) = (k*T(7,3)*(dy^2)+2*k*T(7,1)*(dy^2)+dx*(dx*(h3*Tinf3*dy+k*(T(6,2)+2*T(8,2))+h3*Tinf3*(dy^2)))) / ...
        ((dx^2)*(h3*dy+3*k)+h*dx*(dy^2)+3*k);
    % Lower Right Corner, Node 6
       temps(7,4) = (k*T(7,3)*(dy^2)+(dx^2)*(h3*Tinf3*dy+k*T(8,4))) / ...
           ((dx^2)*(h3*dy+k)+k*(dy^2));
 
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
title("2D Dimensional Steady State Temperatures (Contour)")
xlabel("Length (m)")
ylabel("Height (m)")
c = colorbar;
c.Label.String = "Temperatures °C";

figure(2)
mesh(x,y,FinalTempsSteadyState2D);
title("2D Dimensional Steady State Temperatures (Mesh)")
xlabel("Length (m)")
ylabel("Height (m)")
c = colorbar;
c.Label.String = "Temperatures °C";

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