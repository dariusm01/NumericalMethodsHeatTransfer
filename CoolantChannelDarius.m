%% Adding the functions to the filepath
addpath('ProjectFunctions')

A = zeros(10,5);
b = size(A);
xdim = b(2);
ydim = b(1);

rows = ydim;
cols = xdim;


% Bottom Rectangle
BottomRecX = cm_to_m(0.20); % 0.20cm to m 
BottomRecY = cm_to_m(0.05); % 0.05cm to m 
        
BottomRec_dx = BottomRecX/(xdim-1);
BottomRec_dy = BottomRecY/(ydim-1);

% Left Rectangle
LeftRecX = cm_to_m(0.05); % 0.05cm to m 
LeftRecY = cm_to_m(0.35); % 0.35cm to m 

LeftRec_dx = LeftRecX/(xdim-1);
LeftRec_dy = LeftRecY/(ydim-1);

% Top Rectangle
TopRecX = cm_to_m(0.20); % 0.20cm to m 
TopRecY = cm_to_m(0.1); % 0.1cm to m 

TopRec_dx = TopRecX/(xdim-1);
TopRec_dy = TopRecY/(ydim-1);

% Top Rectangle
for i = 1:2
    for j = 2:cols
        A(i,j) = 10;
    end 
end 

% Left Rectangle
for i = 2:rows-1
    for j = 1
        A(i,j) = 22;
    end 
end 

% Bottom Rectangle
for i = ydim
    for j = 2:cols-1
        A(i,j) = 70;
    end 
end 
      
% Interior 
for i = 3:rows-1
    for j = 2:cols
        A(i,j)= 183;
    end 
end 

% Corners
A(1,1) = 5;
A(1,end) = 5;
A(end,1) = 5;
A(end,end) = 5;