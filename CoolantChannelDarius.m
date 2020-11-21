%% Adding the functions to the filepath
addpath('ProjectFunctions')

%% Setting Up the nodes
A = zeros(50,25);
b = size(A);

%% Rows
TopRecArea = Percent(30); % 30% of the total area
TopRec =  ceil(b(1)*TopRecArea); % 30% of 50

BottomRecArea = Percent(15); % 15% of the total area
BottomRec = ceil(b(1)*BottomRecArea); % 15% of 50 (rounded up)

Remaining = 1-(TopRecArea+BottomRecArea);
CenterRec = floor(b(1)*Remaining); % 55% of 50

CoolantArea = floor((TopRec+CenterRec));

%% Columns
LeftRecAcross = Percent(20); % the left rectangle makes up 20 of the cols
LeftRec = b(2)*LeftRecAcross;
RemainingAcross = 1-LeftRecAcross;
RemainingCols = b(2)*RemainingAcross;

rows = b(1);
cols = b(2);

%% Top Rectangle
for i = 1:TopRec
    for j = LeftRec+1:cols
        A(i,j) = 55;
    end 
end 

%% Left Rectangle
for i = 1:rows
    for j = 1:LeftRec
        A(i,j) = 23;
    end 
end 

%% Center Rec (Coolant Channel)
for i = TopRec:CoolantArea
    for j = LeftRec+1:cols
        A(i,j) = 183;
    end 
end 


%% Bottom Rec 
for i = CoolantArea:(CoolantArea+BottomRec)
    for j = LeftRec+1:cols
        A(i,j) = 97;
    end 
end 


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

%% Iterative Method
iter = 0; % iteration counter
iterLimit = 10000000;

while iter < iterLimit
    %% Outer Corners
    
    % Top Left Corner
    T(1,1) = ((((2*k)/(dx^2)) + ((2*k)/(dy^2)) + ((2*h)/(dy)))^-1) * ((((2*k)/(dx^2))*A(1,2))...
        + (((2*k)/(dy^2))*A(2,1)) + (((2*h)/(dy))*Tinf1) + egen);
    
    % Top Right Corner
    T(1,end) = ((((2*k)/(dx^2)) + ((2*k)/(dy^2)) + ((2*h)/(dy)))^-1) * ((((2*k)/(dx^2))*A(1,end-1))...
        + (((2*k)/(dy^2))*A(2,end)) + (((2*h)/(dy))*Tinf1) + egen);
    
    % Bottom Left Corner
    T(end,1) = ((((2*k)/(dx^2)) + ((2*k)/(dy^2)) + ((2*h_comb)/(dy)))^-1) * ((((2*k)/(dx^2))*A(end,2))...
        + (((2*k)/(dy^2))*A(end-1,1)) + (((2*h_comb)/(dy))*Tinf3) + egen);
    
    % Bottom Right Corner
    T(end,end) = ((((2*k)/(dx^2)) + ((2*k)/(dy^2)) + ((2*h_comb)/(dy)))^-1) * ((((2*k)/(dx^2))*A(end,end-1))...
        + (((2*k)/(dy^2))*A(end-1,end)) + (((2*h_comb)/(dy))*Tinf3) + egen);
    
    %% Inner Corners
    T(TopRec,LeftRec+1) = 855;
    T(TopRec,end) = 855;
    T(CoolantArea-1,LeftRec+1) = 855;
    T(CoolantArea-1,end) = 855;

    
    A = T;
    iter = iter + 1; 
end 


function A = Percent(B)
A = B/100;
end 