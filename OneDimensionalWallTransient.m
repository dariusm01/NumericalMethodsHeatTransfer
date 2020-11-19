%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 2350; % W/m^2k
h2 = 360; % W/m^2k
Tinf1 = KelvintoC(3350); % 3350k to °C
Tinf2 = KelvintoC(295); % 295k to °C
L = cm_to_m(1); % 1cm to m 
n = 10; % number of nodes
dx = L/(n-1); % m [length/(number of nodes - 1)]
rho = 7900; % kg/m^3
cp = Interpolation(300, 200, 477, 402, 295); % J/kg*k
alpha = ThermalDiffusivity(rho, cp, k); % m^2/s
dt = 0.005; % size of steps
tau = meshFourier(alpha, dt,dx); % mesh Fourier number
criteria = 1-2*tau;
timeSteps = 70000; % number of steps

%% Explicit Stability Criterion
if  criteria < 0
    warning("Will not converge, consider decreasing dt")
end 

%% Using the explicit approach
T = zeros(timeSteps, n);
Tinitial =  KelvintoC(295); % 295k to °C
T(1,:) = Tinitial; % setting the first row to the intial temp. These will get updated down the column

for i = 2:timeSteps
    % Boundary 1
    T(i,1) = ((criteria - 2*tau*(h1*dx/k))*T(i-1,1)) + (2*tau*T(i-1,2)) + (2*tau*((h1*dx)/k)*Tinf1) + (tau*egen*dx^2)/k;
    
    % Interior Nodes
    for j = 2:n-1
        T(i,j) = (tau*(T(i-1,j-1)+T(i-1,j+1))) + criteria*T(i-1,j) + (tau*egen*dx^2)/k;
    end 
    
    % Boundary 2
    T(i,n) = ((criteria - 2*tau*(h2*dx/k))*T(i-1,n)) + (2*tau*T(i-1,n-1)) + (2*tau*((h2*dx)/k)*Tinf2) + (tau*egen*dx^2)/k;
end 

FinalTempsTransient1D = T;

%% Plotting Temperature History of boundaries
figure(1)
% All rows, first column
plot(FinalTempsTransient1D(:,1))
grid on
xlabel("Time steps")
ylabel("Temperature °C")
title("Temperature History of Combustion Chamber Boundary")

figure(2)
% All rows, last column
plot(FinalTempsTransient1D(:,n))
grid on
xlabel("Time steps")
ylabel("Temperature °C")
title("Temperature History of Outside Boundary")

function alpha = ThermalDiffusivity(rho, cp, k)
alpha = k/(rho*cp);
end 

function tau = meshFourier(alpha, dt,dx)
tau = (alpha*dt)/(dx^2);
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