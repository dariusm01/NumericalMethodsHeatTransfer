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
iterLimit = 10;

%% Nodes (horizontal & vertical)
dimension = [5 5]; % any # of nodes (x-direction) & nodes (y-direction)

xNodes = dimension(1); % Across
dx = L/(xNodes - 1);

yNodes = dimension(2); % Down
dy = H/(yNodes -1);

T = zeros(xNodes,yNodes); % 'old' temps

temps = zeros(xNodes,yNodes); % 'new' temps

%% 2D Boundaries 
T(1,1:xNodes) = 0;  % Top
T(end,1:xNodes) = 0; % Bottom
T(1:yNodes,1) = 0; % Left
T(1:yNodes,end) = 0; % right




while iter < iterLimit
    
      
            % Upper Left Corner
            
            temps(1,1) = ((h1*dx*(dy^2)*Tinf1) + (k*(dy^2)*T(m+1,n)) + (k*(dx^2)*T(m,n+1)) + (h3*(dx^2)*dy*Tinf3)+...
                ((egen/2)*dx^2*dy^2))/((h1*dx*(dy^2)) + k*dy^2 + k*dx^2 + (h3*(dx^2)*dy));

            % Top
            temps(1,2:xNodes-1) = 1; % Still need eqs 

            % Upper Right Corner
            temps(1,xNodes) = ((h3*dy*(dx^2)*Tinf3) + (k*(dy^2)*T(m-1,n)) + (k*(dx^2)*T(m,n+1)) + (h2*(dy^2)*dx*Tinf2)+...
                ((egen/2)*dx^2*dy^2))/((h3*dy*(dx^2)) + k*dx^2 + k*dy^2 + (h2*(dy^2)*dx));

            % Left Side
            temps(2:yNodes-1,1) = 1;  % Still need eqs 

            % Lower Left Corner
            temps(yNodes,1) = ((h1*dx*(dy^2)*Tinf1) + (k*(dy^2)*T(m+1,n)) + (k*(dx^2)*T(m,n-1)) + (h4*(dx^2)*dy*Tinf4)+...
                ((egen/2)*dx^2*dy^2))/((h1*dx*(dy^2)) + k*dx^2 + k*dy^2 + (h4*(dx^2)*dy));

            % Bottom 
            temps(yNodes,2:xNodes-1) = 1;  % Still need eqs 

            % Lower Right Corner
            temps(xNodes,yNodes) = ((h4*dy*(dx^2)*Tinf4) + (k*(dy^2)*T(m-1,n)) + (k*(dx^2)*T(m,n-1)) + (h2*(dy^2)*dx*Tinf2)+...
                ((egen/2)*dx^2*dy^2))/((h4*dy*(dx^2)) + k*dx^2 + k*dy^2 + (h2*(dy^2)*dx));

            % Right Side
            temps(2:yNodes-1,xNodes) = 1;  % Still need eqs 

        
     % Interior nodes    
    for i = 2:xNodes-1 % rows
        for j = 2:yNodes-1 % colmuns     
        temps(i,j) = ((k*dy^2*(T(i+1,j)+T(i-1,j))) + (k*dx^2*(T(i,j+1)+T(i,j-1))) +...
            (egen*(dx^2)*(dy^2)))/(k*(2*dy^2 +2*dx^2));
        end 
    end 
    
    
    T = temps;
    iter = iter + 1; 
    
end 


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