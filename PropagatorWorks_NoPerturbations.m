clear,clc,close all
% Inputs (your existing values)
mu                  = 1.327e20;
T1                  = datetime(2017,7,1,0,0,0);
T2                  = datetime(2019,12,20);
tolerance           = 1e-5;
Satellite           = bodysort();
Satellite.name      = 'Satellite1';
Satellite.Method    = 'Universal Variable';
Satellite.IttMethod = 'Bisection';
Body_Center         = bodysort();
Body_Center.name    = 'Sun';
Body_Center.mass    = 1.989e30;     % [kg] Sun Mass
Body_Center.radius  = 6.957e8;      % [m]  Sun Radius @ Equator
% Earth:
Body_Start          = bodysort();
Body_Start.name     = 'Earth';
Body_Start.mass     = 5.9722e24;      % [kg] Earth Mass
Body_Start.radius   = 6378e3;         % [m]  Earth Radius @ Equator
% Mars:
Body_Target         = bodysort();
Body_Target.name    = 'Mars';
Body_Target.mass    = 0.64171e24;   % [kg] Mars Mass
Body_Target.radius  = 3396e3;       % [m]  Mars Radius @ Equator

% Obtaining original orbit
[Body_Center, Body_Start, Body_Target, TOF, V1, V2, Satellite, ~] = OriginalOrbit(T1, T2, '421', 'SolarSystem', Body_Center, Body_Target, Body_Start, mu, 3, Satellite, tolerance, 1);


%% Propagation Model 
% Initialise
timer(1)        = juliandate(T1)*(3600*24);
Position(1,:)   = Body_Start.P1;
Velocity(1,:)   = V1;
[~, Orbit_base] = rv2orb(Position(1,:), Velocity(1,:), mu);

% Use the base orbital elements (these stay constant for Keplerian motion)
a = Orbit_base.a;
e = Orbit_base.e;
i = Orbit_base.i;
Omega = Orbit_base.Omega;
w = Orbit_base.w;

% Initialise true anomaly
TA_current(1) = deg2rad(Orbit_base.TA);
TA_increment  = 0.002; % [rad]

target_time = juliandate(T2)*(3600*24);
dt = 60;

% Pre-compute constants
n = sqrt(mu / a^3); % mean motion
E0 = 2 * atan( sqrt((1 - e)/(1 + e)) * tan(TA_current(1)/2) ); % initial E
M0 = E0 - e * sin(E0); % initial mean anomaly
j = 1;

while timer(j) < target_time
    % Advance time
    timer(j+1) = timer(j) + dt;

    % Compute mean anomaly at new time
    M(j) = M0 + n * (timer(j+1) - timer(1));

    % Solve Kepler's equation: M = E - e*sin(E)
    E(j) = solve_kepler_equation(M(j), e, tolerance);

    % Convert E to true anomaly
    TA(j) = 2 * atan2( sqrt(1 + e) * sin(E(j)/2), sqrt(1 - e) * cos(E(j)/2) );
    TA(j) = mod(TA(j), 2*pi); % ensure it's in [0, 2pi]

    % Store and convert to position/velocity
    %TA_current(j+1) = TA(j);
    [Position(j+1,:), Velocity(j+1,:)] = orb2rv(a, e, i, Omega, w, TA(j), mu);

    disp(j)
    j = j + 1;

end

%%
disp(['Total propagation time: ', num2str((timer(end) - timer(1))/86400), ' days']);
disp(['Number of steps: ', num2str(length(timer))]);


%%
figure
plot3(Satellite.PathPosition(:,1), Satellite.PathPosition(:,2), Satellite.PathPosition(:,3), 'k', 'LineWidth', 1);
hold on;
% Plot the perturbed path with a dashed red line
plot3(Position(:,1), Position(:,2), Position(:,3), 'r--', 'LineWidth', 2);
% scatter3(PositionPerturbed(:,1), PositionPerturbed(:,2), PositionPerturbed(:,3), 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');


% % Plot the points with more professional colors and sizes
% plot3(R1(:,1), R1(:,2), R1(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410]); % Blue
% plot3(R2(:,1), R2(:,2), R2(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', [0.8500 0.3250 0.0980]); % Red




%% 







function TA_new = propagate_true_anomaly(TA_current, a, e, mu, dt)
% Propagate true anomaly by time step dt
% Inputs:
%   TA_current - current true anomaly [rad or deg]
%   a - semi-major axis [m]
%   e - eccentricity [-]
%   mu - gravitational parameter [m^3/s^2]
%   dt - time step [s]
%   use_degrees - optional flag, true if angles are in degrees (default: false)
% Output:
%   TA_new - new true anomaly [rad or deg, same units as input]


TA_current_rad = deg2rad(TA_current);


% Step 1: Convert current true anomaly to eccentric anomaly (robust method)
E_current = true_to_eccentric_anomaly(TA_current_rad, e);

% Step 2: Convert eccentric anomaly to mean anomaly
M_current = E_current - e * sin(E_current);

% Step 3: Calculate mean motion
n = sqrt(mu / a^3);  % mean motion [rad/s]

% Step 4: Calculate new mean anomaly
M_new = M_current + n * dt;

% Step 5: Solve Kepler's equation for new eccentric anomaly
E_new = solve_kepler_equation(M_new, e);

% Step 6: Convert new eccentric anomaly to new true anomaly (robust method)
TA_new_rad = eccentric_to_true_anomaly(E_new, e);


TA_new = rad2deg(TA_new_rad);

end

function E = true_to_eccentric_anomaly(TA, e)
% Convert true anomaly to eccentric anomaly using atan2 (handles all quadrants)
% TA: true anomaly [rad]
% e: eccentricity
% E: eccentric anomaly [rad]

% Use atan2 for robust quadrant handling
cos_E = (e + cos(TA)) / (1 + e * cos(TA));
sin_E = sqrt(1 - e^2) * sin(TA) / (1 + e * cos(TA));

E = atan2(sin_E, cos_E);

% Ensure E is in [0, 2π] range
if E < 0
    E = E + 2*pi;
end

end

function TA = eccentric_to_true_anomaly(E, e)
% Convert eccentric anomaly to true anomaly using atan2 (handles all quadrants)
% E: eccentric anomaly [rad]
% e: eccentricity
% TA: true anomaly [rad]

% Use atan2 for robust quadrant handling
cos_TA = (cos(E) - e) / (1 - e * cos(E));
sin_TA = sqrt(1 - e^2) * sin(E) / (1 - e * cos(E));

TA = atan2(sin_TA, cos_TA);

% Ensure TA is in [0, 2π] range
if TA < 0
    TA = TA + 2*pi;
end

end

function E = solve_kepler_equation(M, e, tol)
% Solve Kepler's equation: M = E - e*sin(E) for E
% Uses Newton-Raphson iteration
% Inputs:
%   M - mean anomaly [rad]
%   e - eccentricity [-]
%   tol - tolerance (optional, default 1e-12)

if nargin < 3
    tol = 1e-12;
end

% Normalize M to [0, 2π] range
M = mod(M, 2*pi);

% Initial guess (Danby's method)
if e < 0.8
    E = M + e * sin(M);  % Good for low to moderate eccentricity
else
    E = M + sign(sin(M)) * 0.85 * e;  % Better for high eccentricity
end

% Newton-Raphson iteration
max_iter = 100;
for i = 1:max_iter
    f = E - e * sin(E) - M;           % Function
    df = 1 - e * cos(E);              % Derivative
    
    E_new = E - f / df;               % Newton-Raphson step
    
    if abs(E_new - E) < tol
        E = E_new;
        return;
    end
    
    E = E_new;
end

warning('Kepler equation did not converge after %d iterations', max_iter);
end




function [Body_Center, Body_Start, Body_Target, TOF, V1, V2, Satellite, anomally] = OriginalOrbit(T1, T2, EPH, Center, Body_Center, Body_Target, Body_Start, mu, Dim, Satellite, tolerance, N)
% Function that outputs the original orbit (i.e. solves Lamberts problem
% and gives the orbital elements of the generated orbit
% Inputs:
%   T1          - Time of departure
%   T2          - Time of arrival
%   EPH 
%   Center      - Center of system 
%   Body_Center - Center of transfer
%   Body_Target - Target Body
%   Body_Start  - Start body
%   ICR2_TRANS1 - Transfer Matrix
%   ICR2_TRANS2 - Tranfer Matrix
%   mu          - Standard gravitational parameter
%   Dim         - Dimension
%   Satellite
% Outputs:
%   Body_Center - Data of central body 
%   Body_Start  - Data of starting body
%   Body_Target - Data of target body
%   TOF         - Time of flight
%   V1          - Initial Hyperbolic excess velocity
%   V2          - Final hyperbolic excess velocity
%   Satellite   - Structure including all the spacecraft's parameters
%   anomally    - true anomally

% Invariable Plane Transfer Matrix:
ICR2_TRANS1 = ICRF2IVP(T1,11,EPH,Center);
ICR2_TRANS2 = ICRF2IVP(T2,11,EPH,Center);

% Finding the position of planets of interest 

for body = ['C' 'S' 'T']
    %Grab Data:
    switch body
        case 'C'
            Body = Body_Center;
        case 'S'
            Body = Body_Start;
        case 'T'
            Body = Body_Target;
    end
         
    % Find position and velocity of planets at departure and arrival in ecliptic plane:
    [Body.P1_ecl,Body.V1_ecl] = planetEphemeris(juliandate(T1),Center,Body.name,EPH,'km'); 
    [Body.P2_ecl,Body.V2_ecl] = planetEphemeris(juliandate(T2),Center,Body.name,EPH,'km'); 
         
    % Find position and velocities in invariable plane
    Body.P1= matsolv(Body.P1_ecl,ICR2_TRANS1)*10^3; 
    Body.V1= matsolv(Body.V1_ecl,ICR2_TRANS1)*10^3;
    Body.P2= matsolv(Body.P2_ecl,ICR2_TRANS2)*10^3; 
    Body.V2= matsolv(Body.V2_ecl,ICR2_TRANS2)*10^3;
          
    % Find orbital elements at start and finish
    [Body.E6, Body.Initial] = rv2orb(Body.P1, Body.V1, mu); 
    [Body.E6, Body.Final  ] = rv2orb(Body.P2, Body.V2, mu);
    Body.Avg                = AvgKepElm(Body.Initial, Body.Final);
         
     %Plot Elements:
     [Body.Path1, Body.P11, Body.V11] = orb2pltdta(Dim,Body.Initial,mu);
     [Body.Path2,Body.P21,Body.V21]   = orb2pltdta(Dim,Body.Final,mu);
     %[Body.Path3,~,~]                 = orb2pltdta(Dim,Body.Avg,    Body.radius,mu);
         
       %Save Data:
        switch body
           case 'C'
              Body_Center = Body;
           case 'S'
              Body_Start = Body;
           case 'T'
              Body_Target= Body;
        end
         
end
fprintf('\n\n');
                                                                     
%Radius:
R1 = Body_Start.P1;    %Starting Position
R2 = Body_Target.P2;   %Ending Position
      
%Time of Flight:
[~,~,~,toh,tom,tos] = datevec(between(T1,T2,'time'));          % Time of Flight in hours
TOF                 = 60*(60*toh + tom) + tos;                 % Time of Flight in Seconds
N                   = round((TOF/Body_Target.Initial.TP)-.5);  % Number of Orbits    

% Solving Lambert's Problem:
[V1,V2,~,~,~]   = lambert(mu,R1,R2,TOF,N,Satellite.Method,Satellite.IttMethod,tolerance);

% Finding orbit of satellite:
[~,Satellite.Orbit] = rv2orb(R1,V1,mu);        
         
% Position and velocity along the orbit for all the values of true anomally 
[Satellite.PathPosition,Satellite.PathVelocity, anomally] = orb2pltdta(Dim,Satellite.Orbit,mu);

end