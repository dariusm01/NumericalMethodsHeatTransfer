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
dt = 0.0005; % size of steps
timeSteps = 50; % number of steps

%% Nodes (horizontal & vertical)
dimension = [5 5]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

xNodes = dimension(2); % Across
dx = L/(xNodes - 1);

yNodes = dimension(1); % Down
dy = H/(yNodes - 1);

rows = yNodes;
cols = xNodes;


%% Using the explicit approach
n = rows*cols; % total number of nodes
T = zeros(timeSteps, n);
Tinitial =  KelvintoC(295); % 295k to °C
T(1,:) = Tinitial; % setting the first row to the intial temp. These will get updated down the column


TwoDNodes = NodeSystem(rows, cols);
temps = zeros(size(TwoDNodes));
Ts = 0.5*ones(size(TwoDNodes)); % very generic

%% Referring to Node system
% Top = TwoDNodes(1,2:yNodes-1)
% Bottom = TwoDNodes(end,2:yNodes-1)
% Left side = TwoDNodes(2:xNodes-1,1)
% Right Side = TwoDNodes(2:xNodes-1,end)
for k = 2:timeSteps
    %% Populating the corners
    % Upper Left corner 
    temps(1,1) = Ts(1,1) + ((-Ts(1,1)*(dx^2 + dx^2 + (((h3*(dx^2)*dy) + (h1*(dy^2)*dx))/k))) + (dy^2*(Ts(1,2)))+...
        (dx^2*(Ts(2,1))) + (((h3*(dx^2)*dy)/k)*Tinf3) + (((h1*(dy^2)*dx)/k)*Tinf1) + ((egen*(dx^2)*(dy^2))/(2*k)))*...
        ((2*alpha*dt)/((dx^2)*(dy^2)));  

    % Bottom Left corner 
    temps(end,1) = Ts(end,1) + ((-Ts(end,1)*(dx^2 + dx^2 + (((h4*(dx^2)*dy) + (h1*(dy^2)*dx))/k))) + (dy^2*(Ts(end,2)))+...
        (dx^2*(Ts(end-1,1))) + (((h4*(dx^2)*dy)/k)*Tinf4) + (((h1*(dy^2)*dx)/k)*Tinf1) + ((egen*(dx^2)*(dy^2))/(2*k)))*...
        ((2*alpha*dt)/((dx^2)*(dy^2))); 

    % Upper Right corner
    temps(1,end) = Ts(1,end) + ((-Ts(1,end)*(dx^2 + dx^2 + (((h3*(dx^2)*dy) + (h2*(dy^2)*dx))/k))) + (dy^2*(Ts(1,end-1)))+...
        (dx^2*(Ts(2,end))) + (((h3*(dx^2)*dy)/k)*Tinf3) + (((h2*(dy^2)*dx)/k)*Tinf2) + ((egen*(dx^2)*(dy^2))/(2*k)))*...
        ((2*alpha*dt)/((dx^2)*(dy^2)));

    % Bottom Right corner 
    temps(end,end) = Ts(end,end) + ((-Ts(end,end)*(dx^2 + dx^2 + (((h4*(dx^2)*dy) + (h2*(dy^2)*dx))/k))) + (dy^2*(Ts(end,end-1)))+...
        (dx^2*(Ts(end-1,end))) + (((h4*(dx^2)*dy)/k)*Tinf4) + (((h2*(dy^2)*dx)/k)*Tinf2) + ((egen*(dx^2)*(dy^2))/(2*k)))*...
        ((2*alpha*dt)/((dx^2)*(dy^2)));  


%     %% Populating the top and bottom
%     for i=2:rows-1
%         for j = 2:cols-1
% 
%             % Top
%             if TwoDNodes(i,j) > TwoDNodes(1) && TwoDNodes(i,j) < TwoDNodes(1,end)
%                 temps(i,j) = ((((alpha*dt)/(dx^2))/dy)*((h3*dx*dy*(Tinf2-Ts(1,j))/k)) + (dx^2*(T(2,j)-T(1,j))/dy)...
%             + dy*(T(1,j-1)+T(1,j+1)-2*T(1,j))/2 + egen*dx^2*dy/2*k) + T(1,j); % generic for now, insert eq     
% 
%             % Bottom    
%             elseif TwoDNodes(i,j) > TwoDNodes(end,1) && TwoDNodes(i,j) < TwoDNodes(end,end)
%                 temps(i,j) = ((((alpha*dt)/(dx^2))/dy)*((h4*dx*dy*(Tinf2-Ts(end,j))/k)) + (dx^2*(T(end-1,j)-T(end,j))/dy)...
%             + dy*(T(end,j-1)+T(end,j+1)-2*T(end,j))/2 + egen*dx^2*dy/2*k) + T(end,j); % generic for now, insert eq           
%             end 
%         end 
%     end 

    %% Populating the sides
    for i = 2:rows-1

        % Left side
        temps(i,1) = ((((alpha*dt)/(dx^2))/dy)*((h1*dx^2*(Tinf1-Ts(i,1))/k)) + (dy*(T(i,2)-T(i,1)))...
            + dx^2*(T(i-1,1)+T(i+1,1)-2*T(i,1))/2*dy + egen*dx^2*dy/2*k) + T(i,1); % generic for now, insert eq     

        % Right side
        temps(i,end) = ((((alpha*dt)/(dx^2))/dy)*((h2*dx^2*(Tinf2-Ts(i,end))/k)) + (dy*(T(i,end-1)-T(i,end)))...
            + dx^2*(T(i-1,end)+T(i+1,end)-2*T(i,end))/2*dy + egen*dx^2*dy/2*k) + T(i,end); % generic for now, insert eq    
       
        %Top
        temps(i,j) = ((((alpha*dt)/(dx^2))/dy)*((h3*dx*dy*(Tinf2-Ts(1,j))/k)) + (dx^2*(T(2,j)-T(1,j))/dy)...
             + dy*(T(1,j-1)+T(1,j+1)-2*T(1,j))/2 + egen*dx^2*dy/2*k) + T(1,j); % generic for now, insert eq 
         
         %Bottom
         temps(i,j) = ((((alpha*dt)/(dx^2))/dy)*((h4*dx*dy*(Tinf2-Ts(end,j))/k)) + (dx^2*(T(end-1,j)-T(end,j))/dy)...
             + dy*(T(end,j-1)+T(end,j+1)-2*T(end,j))/2 + egen*dx^2*dy/2*k) + T(end,j); % generic for now, insert eq
    end 

    %% Interior Nodes
    for i = 2:rows-1
        for j = 2:cols-1

            temps(i,j) = (((alpha*dt)/((dx^2)*(dy^2))) * ((-Ts(i,j)*(2*dy^2 + 2*dx^2)) + (dy^2*(Ts(i,j-1)+Ts(i,j+1)))+... 
                (dx^2*(Ts(i-1,j)+Ts(i+1,j))) + (egen*(dx^2 * dy^2)/k))) + Ts(i,j);
        end 
    end 

    Ts = temps;

B = reshape(temps,[1,n]);

T(k,:) = B;
end 

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
