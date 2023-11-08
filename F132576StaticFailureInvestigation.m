%% Variable Declaration
% F132576
%A = 576 [Default 0];
A = 576;
%B = 132*576*10 [Default 100*1000*10];
B = 132*576*10;

% Forces_ABC[Y comp, Z comp]
Forces_A        = [2500+A, 2000]
Forces_B        = [2500+A, 2000]
Forces_C        = [-5000-(2*A),-4000]

% Other Material
n               = 2                       % Arbitarily set.
Sy              = 380+B/(1e6)             % As given in Coursework
tau_allowable   = (Sy/(2*n))              % Maximum shear theory design equation
tau_dash        = Sy/n                    % Yield design equation (Von-Mises)

% Distances (mm)
CA              = 105;                    % Distance from C to A
CB              = CA;                     % Distance from C to B, isn't really used.
DB              = 37;                     % Distance from D to B
Ax_AC           = 55;                     % Axis to axis

% Other params
I               = pi/64;                  % in terms of d^4.
J               = 2*I;

% Global Torque
SystemTorque    = abs(Forces_C(2)*Ax_AC)

%% Investigation
% Investigating the bending moment from C to A
Moment_CZ        = Forces_A(2)*abs(CA)    % Moment at C due to A in the Z Direction
Moment_CY        = Forces_A(1)*abs(CA)    % Moment at C due to A in the Y Direction

sigma_3Coeff     = (Moment_CZ*(0.5))/(I)  % Over d^3
sigma_1Coeff     = (Moment_CY*(0.5))/(I)  % Over d^3

tau_maxCoeff1    = abs((sigma_1Coeff - sigma_3Coeff)/2) % Tau_Max Equation (Over d^3)

d1              = nthroot(tau_maxCoeff1/tau_allowable,3)

%d1 = (2*nthroot(4,3)*nthroot((Moment_C),3))/ ...
    %(nthroot(pi*sigma_b,3))
  

% Investigating the torsion at D due to B

Moment_D        = [Forces_B(1)*abs(DB);             % Y
                   Forces_B(2)*abs(DB)]             % Z
sigma_b         = [(Moment_D(1)*(0.5))/(I);         % Y
                   (Moment_D(2)*(0.5))/(I)]         % Z (Over d^3 in calculations...)

tau_torsion     = SystemTorque*0.5/J                % (Over d^3)

tau_maxCoeff2   = sqrt( ...
                  ( ((abs(sigma_b(1)-sigma_b(2)))/2) ^2) + ...
                  (tau_torsion)^2)                  
                                                    % Big fat equation, radius for tau-max
                                                    % Absolutes are used
                                                    % here as to not break
                                                    % the shy equation.
                                                    
%disp("Tau_max 2: "+ tau_maxCoeff2+"/d^3")  %DEBUGGING!

d2              = nthroot(tau_maxCoeff2/tau_allowable,3)

%% D due to B: Re-calculation using von-mises equations.
sigma1          = sigma_b(1)        %
sigma3          = sigma_b(2)        % 
sigma_vonmis    = (sqrt(sigma1^2 + sigma3^2 - (sigma1*sigma3))) % Over d^3
d3              = nthroot(sigma_vonmis/tau_dash,3) % Solve for d
% Doesn't make much sense as torsion isn't considered?

%% Finalisation
finalD          = max([d3,d2,d1]);
disp("Checking both C & D and their respective diameters,")
disp("The minimum diameter for the shaft should be "+finalD)
disp("However, in the context of machine elements, it should be preferred to accept a rounded or standardised number, therefore d may be set to "+round(finalD)+ " or whichever closest size following a metrologic standard.")

%% Mohr's Circle (Optional):
% Define principal stresses
sigma1 = sigma_3Coeff; % Replace with your value
sigma3 = sigma_1Coeff;  % Replace with your value

% Calculate the center and radius of Mohr's circle
center = (sigma1 + sigma3) / 2;
radius = (sigma1 - sigma3) / 2;

% Calculate the normal and shear stresses at an angle theta (in radians)
theta = linspace(0, pi, 100); % You can adjust the range and number of points
sigma_n = center + radius * cos(2 * theta);
tau = radius * sin(2 * theta);

% Plot Mohr's circle
figure;
plot(sigma_n, tau, 'b', 'LineWidth', 2);
title('Mohr''s Circle');
xlabel('\sigma_n (Normal Stress)');
ylabel('\tau (Shear Stress)');
grid on;

% Add principal stresses as points on the circle
hold on;
plot(sigma1, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(sigma3, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
legend('\sigma_n - \tau Plot', '\sigma_1', '\sigma_3');

% Display the results
fprintf('Center of Mohr''s circle: %.2f\n', center);
fprintf('Radius of Mohr''s circle: %.2f\n', radius);

% Optional: You can also add a grid for clarity
grid on;

