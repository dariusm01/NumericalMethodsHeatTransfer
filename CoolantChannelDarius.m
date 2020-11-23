%% Adding the functions to the filepath
addpath('ProjectFunctions')

%% Nodes (horizontal & vertical)
dimension = [50 25]; % any # of nodes (x-direction) & nodes (y-direction)
% similar to a coordinate (x,y)

%% Setting Up the nodes
A = zeros(dimension(1) ,dimension(2));

%% Rows
TopRecArea = Percent(30); % 30% of the total area
TopRec =  ceil(dimension(1)*TopRecArea); % 30% of 50

BottomRecArea = Percent(15); % 15% of the total area
BottomRec = ceil(dimension(1)*BottomRecArea); % 15% of 50 (rounded up)

Remaining = 1-(TopRecArea+BottomRecArea);
CenterRec = floor(dimension(1)*Remaining); % 55% of 50

CoolantArea = floor((TopRec+CenterRec));

%% Columns
LeftRecAcross = Percent(20); % the left rectangle makes up 20 of the cols
LeftRec = dimension(2)*LeftRecAcross;
RemainingAcross = 1-LeftRecAcross;
RemainingCols = dimension(2)*RemainingAcross;

rows = dimension(1);
cols = dimension(2);


T = zeros(size(A));

%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h = 360; % W/m^2k Air
h_cool = 8800; % W/m^2k Coolant Channel
h_comb = 2350; % W/m^2k Combustion Chamber
Tinf1 = KelvintoC(295); % 295k to °C (Air)
Tinf2 = KelvintoC(183); % 183k to °C (Coolant Channel)
Tinf3 = KelvintoC(3350); % 3350k to °C (Combustion Chamber)
L = cm_to_m(0.25); % 1cm to m 
H = cm_to_m(0.35); % 30cm to m 

TopRecLength = cm_to_m(0.20);
TopRecHeight = cm_to_m(0.10);

CenterRecLength = TopRecLength;
CenterRecHeight = 2*TopRecHeight;

BottomRecLength = TopRecLength;
BottemRecHeight = 0.5*TopRecHeight;

LeftRecLength = BottemRecHeight;
LeftRecHeight = TopRecHeight + CenterRecHeight + BottemRecHeight;

%% Different dy's and dx's

dx_topRec = TopRecLength/(RemainingCols-1);
dy_topRec = TopRecHeight/(TopRec-1);

dx_centerRec = CenterRecLength/(RemainingCols-1);
dy_centerRec = CenterRecHeight/(CenterRec-1);

dx_bottomRec = BottomRecLength/(RemainingCols-1);
dy_bottomRec = BottemRecHeight/(BottomRec-1);

dx_leftRec = LeftRecLength/((dimension(2)-RemainingCols)-1);
dy_leftRec = LeftRecHeight/(dimension(1)-1);

%% Iterative Method
iter = 0; % iteration counter
iterLimit = 10000;

while iter < iterLimit
    %% Outer Corners
    
    % Top Left Corner
    T(1,1) = ((((2*k)/(dx_leftRec^2)) + ((2*k)/(dy_leftRec^2)) + ((2*h)/(dy_leftRec)))^-1) * ((((2*k)/(dx_leftRec^2))*A(1,2))...
        + (((2*k)/(dy_leftRec^2))*A(2,1)) + (((2*h)/(dy_leftRec))*Tinf1) + egen);
    
    % Top Right Corner
    T(1,end) = ((((2*k)/(dx_topRec^2)) + ((2*k)/(dy_topRec^2)) + ((2*h)/(dy_topRec)))^-1) * ((((2*k)/(dx_topRec^2))*A(1,end-1))...
        + (((2*k)/(dy_topRec^2))*A(2,end)) + (((2*h)/(dy_topRec))*Tinf1) + egen);
    
    % Bottom Left Corner
    T(end,1) = ((((2*k)/(dx_leftRec^2)) + ((2*k)/(dy_leftRec^2)) + ((2*h_comb)/(dy_leftRec)))^-1) * ((((2*k)/(dx_leftRec^2))*A(end,2))...
        + (((2*k)/(dy_leftRec^2))*A(end-1,1)) + (((2*h_comb)/(dy_leftRec))*Tinf3) + egen);
    
    % Bottom Right Corner
    T(end,end) = ((((2*k)/(dx_bottomRec^2)) + ((2*k)/(dy_bottomRec^2)) + ((2*h_comb)/(dy_bottomRec)))^-1) * ((((2*k)/(dx_bottomRec^2))*A(end,end-1))...
        + (((2*k)/(dy_bottomRec^2))*A(end-1,end)) + (((2*h_comb)/(dy_bottomRec))*Tinf3) + egen);
    
    %% Inner Corners

    % Inner top left corner
    T(TopRec,LeftRec+1) = (((4*k)/(dy_topRec^2)) + ((2*h_cool)/dy_topRec) + ((2*k)/(dx_topRec^2)))^-1 * (((4*k)/(dy_topRec^2)*A(TopRec,LeftRec)) +...
        (((2*h_cool)/dy_topRec)*Tinf2) + (k/(dx_topRec^2))*(A(TopRec,LeftRec) + A(TopRec,LeftRec+2)) + (3*egen));
    
    % Inner top right corner
    T(TopRec,end) = (((2*k)/(dy_topRec^2)) + ((2*h_cool)/dy_topRec) + ((2*k)/(dx_topRec^2)))^-1 * (((2*k)/(dy_topRec^2)*A(TopRec,end-1)) +...
        (((2*h_cool)/dy_topRec)*Tinf2) + ((2*k)/(dx_topRec^2))*(A(TopRec,end-1)) + (egen));
    
    % Inner bottom left corner
    T(CoolantArea-1,LeftRec+1) = (((4*k)/(dy_bottomRec^2)) + ((2*h_cool)/dy_bottomRec) + ((2*k)/(dx_bottomRec^2)))^-1 * (((4*k)/(dy_bottomRec^2)*A(CoolantArea,LeftRec+1)) +...
        (((2*h_cool)/dy_bottomRec)*Tinf2) + (k/(dx_bottomRec^2))*(A(CoolantArea-1,LeftRec) + A(CoolantArea+1,LeftRec+2)) + (3*egen));
    
    % Inner bottom right corner
    T(CoolantArea-1,end) = (((2*k)/(dy_bottomRec^2)) + ((2*h_cool)/dy_bottomRec) + ((2*k)/(dy_bottomRec^2)))^-1 * (((2*k)/(dy_bottomRec^2)*A(CoolantArea,end)) +...
        (((2*h_cool)/dy_bottomRec)*Tinf2) + ((2*k)/(dx_bottomRec^2))*(A(CoolantArea-1,end-1)) + (egen));

    
    % Left side
    for m = 2:rows-1
        for n = 2:TopRec-1
        T(m,1) = ((k*dx_leftRec^2)*(A(m-1,n) + A(m+1,n)) + (2*k*((dy_leftRec^2)*A(m,n+1))) + (egen*(dx_leftRec^2*dy_leftRec^2)))/(2*k*(dx_leftRec^2 + dy_leftRec^2));
        
    % Left side Interior
       T(m,n) = ((k*dy_leftRec^2*(A(m,n+1)+A(m,n-1))) + (k*dx_leftRec^2*(A(m+1,n)+ A(m-1,n))) +...
        (egen*(dx_leftRec^2)*(dy_leftRec^2)))/(k*(2*dy_leftRec^2 +2*dx_leftRec^2));
        end 
    end 
    
    %% Right side 
    
    for m = 2:TopRec-1
        for n = LeftRec+1:cols-1
    % Right side top rec
    T(m,end) = ((k*(dx_topRec^2)*(A(m-1,n) + A(m+1,n))) +...
            (2*k*(dy_topRec^2)*A(m,n+1)) + (egen*(dx_topRec^2)*(dy_topRec^2)))/((2*k*(dx_topRec^2)) +...
            (2*k*(dy_topRec^2)));
        
     % Top rec interior  
     T(m,n) = ((k*dy_topRec^2*(A(m,n+1)+A(m,n-1))) + (k*dx_topRec^2*(A(m+1,n)+A(m-1,n))) +...
        (egen*(dx_topRec^2)*(dy_topRec^2)))/(k*(2*dy_topRec^2 +2*dx_topRec^2));   
        end 
    end
        
    for m = TopRec:CoolantArea-1
        for n = LeftRec+1:cols-1   
    % Right side center rec    
    T(m,end) = ((2*h_cool*(dy_centerRec^2)*dx_centerRec*Tinf1) + (k*(dx_centerRec^2)*(A(m-1,n) + A(m+1,n))) +...
        (2*k*(dy_centerRec^2)*A(m,n+1)) + (egen*(dx_centerRec^2)*(dy_centerRec^2)))/((2*h_cool*(dy_centerRec^2)*dx_centerRec) +...
        (2*k*(dx_centerRec^2)) + (2*k*(dy_centerRec^2)));  
    
    % Center rec interior
    T(m,n) = ((k*dy_centerRec^2*(A(m,n+1)+A(m,n-1))) + (k*dx_centerRec^2*(A(m+1,n)+A(m-1,n))) +...
        (egen*(dx_centerRec^2)*(dy_centerRec^2)))/(k*(2*dy_centerRec^2 +2*dx_centerRec^2));
        end 
    end 
        
    for m = CoolantArea:(CoolantArea+BottomRec)-1
        for n = LeftRec+1:cols-1
    % Right side bottom rec
    T(m,end) = ((k*(dx_bottomRec^2)*(A(m-1,n) + A(m+1,n))) +...
            (2*k*(dy_bottomRec^2)*A(m,n+1)) + (egen*(dx_bottomRec^2)*(dy_bottomRec^2)))/((2*k*(dx_bottomRec^2)) +...
            (2*k*(dy_bottomRec^2)));
    
    % Bottom rec interior
    T(m,n) = ((k*dy_bottomRec^2*(A(m,n+1)+A(m,n-1))) + (k*dx_bottomRec^2*(A(m+1,n)+A(m-1,n))) +...
        (egen*(dx_bottomRec^2)*(dy_bottomRec^2)))/(k*(2*dy_bottomRec^2 +2*dx_bottomRec^2));
        end 
    end 

    %% Top
    
    % Left Rectangle
    for m = 1
        for n = 2:LeftRec-1
            T(1,n) = ((2*h*(dx_leftRec^2)*dy_leftRec*Tinf1) + (k*(dy_leftRec^2)*(A(m,n-1) + A(m,n+1))) +...
                (2*k*(dx_leftRec^2)*A(m+1,n)) + (egen*(dx_leftRec^2)*(dy_leftRec^2)))/((2*h*(dx_leftRec^2)*dy_leftRec) + (2*k*(dy_leftRec^2)) +...
                (2*k*(dx_leftRec^2)));
        end 
    end 
    
    % Top Rectangle
    for m = 1
        for n = LeftRec:cols-1
            T(1,n) = ((2*h*(dx_topRec^2)*dy_topRec*Tinf1) + (k*(dy_topRec^2)*(A(m,n-1) + A(m,n+1))) +...
                (2*k*(dx_topRec^2)*A(m+1,n)) + (egen*(dx_topRec^2)*(dy_topRec^2)))/((2*h*(dx_topRec^2)*dy_topRec) + (2*k*(dy_topRec^2)) +...
                (2*k*(dx_topRec^2)));
        end 
    end 

    %% Bottom
    
    % Left Rectangle
    for m = rows
        for n = 2:LeftRec-1
        T(rows,n) = ((2*h_comb*(dx_leftRec^2)*dy_leftRec*Tinf3) + (k*(dy_leftRec^2)*(A(m,n-1) + A(m,n+1))) +...
            (2*k*(dx_leftRec^2)*A(m-1,n)) + (egen*(dx_leftRec^2)*(dy_leftRec^2)))/((2*h_comb*(dx_leftRec^2)*dy_leftRec) + (2*k*(dx_leftRec^2)) +...
            (2*k*(dx_leftRec^2)));
        end 
    end 
    
    % Bottom Rec 
    for m = rows
        for n = LeftRec:cols-1
            T(rows,n) = ((2*h_comb*(dx_bottomRec^2)*dy_bottomRec*Tinf3) + (k*(dy_bottomRec^2)*(A(m,n-1) + A(m,n+1))) +...
            (2*k*(dx_bottomRec^2)*A(m-1,n)) + (egen*(dx_bottomRec^2)*(dy_bottomRec^2)))/((2*h_comb*(dx_bottomRec^2)*dy_bottomRec) + (2*k*(dy_bottomRec^2)) +...
            (2*k*(dx_bottomRec^2)));
        end 
    end 
        
    A = T;
    iter = iter + 1; 
end 

FinalTempsCoolantChannel = T;

%% Plotting
figure(1)
contourf(FinalTempsCoolantChannel);
title("Coolant Channel Steady State Temperatures (Contour)")
xlabel("Length (m)")
ylabel("Height (m)")
c = colorbar;
c.Label.String = "Temperatures °C";

figure(2)
mesh(FinalTempsCoolantChannel);
title("Coolant Channel Steady State Temperatures (Mesh)")
xlabel("Length (m)")
ylabel("Height (m)")
c = colorbar;
c.Label.String = "Temperatures °C";