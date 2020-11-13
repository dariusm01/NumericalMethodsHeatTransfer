%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 2350; % W/m^2k
h2 = 360; % W/m^2k
Tinf1 = KelvintoC(3350); % 3350k to °C
Tinf2 = KelvintoC(295); % 295k to °C
L = 1e-2; % m 
n = 100; % number of nodes
dx = L/(n-1); % m [length1/(number of nodes - 1)]
iter = 0; % iteration counter
iterLimit = 1000000;
x = 0:dx:L;


%% Using Direct Method

A = zeros(n,n);
b = zeros(n,1);

%BC Left
A(1,1) = -(h1*dx +k);
A(1,2) = k;
b(1,1) = -Tinf1*h1*dx - (egen*dx^2)/2;

% Interior Points
for i = 2:n-1
    A(i,i-1) = 1;
    A(i,i) = -2;
    A(i,i+1) = 1;
    b(i,1) = -(egen/k)*dx^2;
end

%BC Right
A(n,n-1) = k;
A(n,n) = -(h2*dx+k);
b(n,1) = -Tinf2*h2*dx - (egen*dx^2)/2;

T2 = A\b;


function T = KelvintoC(x)
T = x-273.15;
end 

function x = Interpolation(y2, y1, x2, x1, YourVal)

m = (y2-y1)/(x2-x1);

x = ((YourVal - y1)/m) + x1;
end 