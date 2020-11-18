%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 2350; % W/m^2k (left)
h2 = 360; % W/m^2k (right)
h3 = 5; % W/m^2k (top) 
h4 = 45; % W/m^2k (bottom)
Tinf1 = KelvintoC(3350); % 3350k to °C
Tinf2 = KelvintoC(295); % 295k to °C
L = cm_to_m(1); % 1cm to m 
H = cm_to_m(30); % 30cm to m 
iter = 0; % iteration counter
iterLimit = 1000000;

%% Nodes (horizontal & vertical)
dimension = [5 5]; % any # of nodes (x-direction) & nodes (y-direction)

xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes -1);

T = zeros(xNodes,yNodes);

%% 2D Boundaries 
T(1,1:xNodes) = 0;  % Top
T(end,1:xNodes) = 0; % Bottom
T(1:yNodes,1) = 0; % Left
T(1:yNodes,end) = 0; % right

%% Corners 
T(1,1) = (h1*Tinf1/dy + k*T(m+1,n)/dx^2 + h3*Tinf2/dx + k*T(m,n+1)/dy^2 + egen/2)/(h1/dy + k/dx^2 + h3/dx + k/dy^2); % Upper Left Corner
%Using a 5x5 Box 
T(end,1) = (k*T(m-1,n)/dx^2 + Tinf2*(h2/dy + h3/dx) + k*T(m,n+1)/dy^2 + egen/2)/(k/dx^2 + h2/dy + h3/dx + k/dy^2); % Upper Right Corner
T(1,end) = (h1*Tinf1/dy + k*T(m+1,n)/dx^2 + h4*Tinf2/dx + k*T(m,n-1)/dy^2 + egen/2)/(h1/dy + k/dx^2 + h4/dx + k/dy^2); % Lower Left Corner
T(end,end) = (k*T(m-1,n)/dx^2 + Tinf2*(h2/dy +h4/dx) + k*T(m,n-1)/dy^2 + egen/2)/(k/dx^2 + h2/dy + h4/dx +k/dy^2); % Lower Right Corner




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