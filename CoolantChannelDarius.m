%% Adding the functions to the filepath
addpath('ProjectFunctions')

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


%% Outer Corners
A(1,1) = 987;
A(1,end) = 987;
A(end,1) = 987;
A(end,end) = 987;

%% Inner Corners
A(TopRec,LeftRec+1) = 855;
A(TopRec,end) = 855;
A(CoolantArea-1,LeftRec+1) = 855;
A(CoolantArea-1,end) = 855;

function A = Percent(B)
A = B/100;
end 