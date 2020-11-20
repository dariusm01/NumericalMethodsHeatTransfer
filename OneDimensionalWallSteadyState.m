%% Defining paramters
k = Interpolation(300, 200, 14.9, 12.6, 295); % W/mk
egen = 0; % W/m^3
h1 = 2350; % W/m^2k
h2 = 360; % W/m^2k
Tinf1 = KelvintoC(3350); % 3350k to °C
Tinf2 = KelvintoC(295); % 295k to °C
L = cm_to_m(1); % 1cm to m 
n = 100; % number of nodes
dx = L/(n-1); % m [length1/(number of nodes - 1)]
iter = 0; % iteration counter
iterLimit = 1000000;
x = 0:dx:L;

%% Using the iterative approach
T = zeros(1,n); % 'old' temps

temps = zeros(1,n); % 'New' temps

while iter < iterLimit
    
    % Boundary 1
    temps(1) = ((Tinf1*h1*dx + k*T(2) + egen*((dx^2)/2))/(h1*dx + k)); 
    
    % Interior Nodes
    for i = 2:n-1
        temps(i) = ((egen/k) + (T(i-1)/(dx^2)) + (T(i+1)/(dx^2)))*((dx^2)/2);
    end
    
    % Boundary 2
    temps(n) = ((Tinf2*h2*dx + k*T(n-1) + egen*((dx^2)/2))/(h2*dx + k)); 
    
    % redefining for next iteration
    T = temps;
    iter = iter + 1; 
end 

FinalTempsSteadyState1D = temps;

%% Plotting
figure(1)
plot(x,FinalTempsSteadyState1D)
grid on
title("1 Dimensional Steady State Temperatures")
xlabel("Length (m)")
ylabel("Temperatures °C")

