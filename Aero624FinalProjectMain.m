% Aero 624 Final Project
% The overall goal of this project is too code a simulation that relicates
% the J2 effect that earth has on a sattelite's orbital elements in MATLAB,
% using Gauss's Planetary Equations and then produce the same result using
% the Unreal Engine. UE4 is the platform that our research group, 
% SpaceCRAFT, advised by Dr. Chamitoff, uses to code our simulations.

% I am starting with MATLAB because I have become very familiar with the
% syntax and nuances of using the built in numerical solvers like ode45.
% All of the simulations within the unreal engine are coded in C++ and
% don't have such luxuries, so a large task of the project will be
% learning how to solve differential equations in C++.

% With that said, there are many possible applications to including this
% functionality to SpaceCRAFT and it's why I've chosen this particular
% task. If a user wants to simulate how a particular orbit configuration
% will affect their satellite, they can simply "plug in" the functionality
% of the code and the effect that J2 has can be analyzed in the scope of
% their mission.

% The format of the code will be set up to where the user can specify the
% mission and satellite parameters, and then visualize what the orbit of
% the satellite will look like over time. The initial tests will use orbits
% with a high enough altitude that drag has negligible effects.

% Future iterations of the code will include atmospheric drag effects to
% increase the fidelity of the simulation. I also plan on incorporating the
% third body problem so that users can begin to plan interplanetary
% trajectories. This is a very crucial aspect to the what the
% functionalities should be like iin SpaceCRAFT and why I think it's a good
% option to choose this as my final Project. Do note, that what I turn in
% will not be the final product, work will continue over the summer and
% into the fall with hopes that a robust system can be included in a
% software release at the end of the Fall semester or in the spring.


clear; clc; close all;

global mu Re J2

mu = 398601.2;
Re = 6371; %Km
J2 = 1.081874E-3; 

% Setting up the user specified orbit parameters
% Radius at Perigee: Recommended value > Re + 1000 km

%%%%  Default parameters are as follows: %%%%
r_p0      = Re + 300;
r_a0      = Re + 3000;
Omega0    = deg2rad(45);
inc0      = deg2rad(50);
argp0     = deg2rad(30);
TA0       = deg2rad(40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prompt user to select default values or custom
DefaultYN = input('Would you like to use default values (y/n):', 's');


if DefaultYN == 'y'
    %%%%  Default parameters are as follows: %%%%
    h_p0      = 1000;
    h_a0      = 3000;
    Omega_deg = 30;
    inc_deg   = 50;
    argp_deg  = 20;
    TA_deg    = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
    
    h_p0 = input('Initial altitude at Perigee (a value of 1000+ is recommended):  ');
    r_p0 = h_p0 + Re;

    % Radius at Apogee: Recommended value > Re + 1000 km 
    h_a0 = input('Initial altitude at Apogee: ');

    while h_a0 < h_p0
        h_a0 = input('Please select value greater than Radius at Perigee: ');
    end 
    r_a0 = h_a0 + Re;

    % Prompt user to input initial RA
    Omega_deg = input('Please input Right Ascension of the Ascending Node parameter (Degrees between 0 and 360): ');

    while Omega_deg > 360 || Omega_deg < 0
        input('Please input a value between 0 and 2*pi: ');
    end
    Omega0 = deg2rad(Omega_deg);

    % Prompt the user for inclination
    inc_deg = input('Please input a value for inclination (degrees. 0 to 90): ');

    while inc_deg > 90 || inc_deg < 0
        input('Please input a value between 0 and 90: ');
    end
    inc0 = deg2rad(Omega_deg);

    % Prompt user for argument of perigee
    argp_deg = input('Please input Argument of Perigee parameter (Degrees between 0 and 360): ');

    while argp_deg > 360 || Omega_deg < 0
        input('Please input a value between 0 and 2*pi: ');
    end
    argp0 = deg2rad(argp_deg);

    %Prompt user for true anamoly
    TA_deg = input('Please input True Anomaly parameter (Degrees between 0 and 360): ');
    while TA_deg > 360 || TA_deg < 0
        input('Please input a value between 0 and 2*pi: ');
    end
    TA0 = deg2rad(TA_deg);
    
end

% Calculated Keplerian orbital elements
e0 = (r_a0 - r_p0)/(r_p0 + r_a0);
h0 = sqrt(r_p0*mu*(1 + e0));
a = (r_a0 + r_p0)/2;
n = sqrt(a^3/mu);
T = 2*pi*n;

% Storing orbital elements
Kep = [h0 e0 Omega0 inc0 argp0 TA0]';

% Prompt user for number of orbits
NumOrbits = input('Please input number of orbits: ');
NumPoints = 500;
tspan = linspace(0, NumOrbits*T, NumPoints*NumOrbits);

eps1 = 1E-12;
options = ['reltol', eps1, 'abstol', eps1];

[tout, dkepdt] = ode45('ODEJ2Gauss', tspan, Kep, options);

h = dkepdt(:,1);
e = dkepdt(:,2);
Omega = dkepdt(:,3);
inc = dkepdt(:,4);
argp = dkepdt(:,5);
TA = dkepdt(:,6);

figure(1)
subplot(3,2,1)
plot(tspan/T, h)
title('Change in Angular Momentum Vector (h) over time')

subplot(3,2,2)
plot(tspan/T, e)
title('Change in Eccentricity (e) over time')

subplot(3,2,3)
plot(tspan/T, Omega)
title('Change in Right Ascension of the Ascending Node (\Omega) over time')

subplot(3,2,4)
plot(tspan/T, inc)
title('Change in Inclination (i) over time')

subplot(3,2,5)
plot(tspan/T, argp)
title('Change in Argument of Perigee (\omega) over time')

subplot(3,2,6)
plot(tspan/T, TA)
title('Change in True Anomaly (\Theta) over time')

% Converting the changing orbital elements to position and velocity at each
% time step

rbarnew = [];
vbarnew = [];
for i = 1:length(tspan)
    [rbarnew(i,:), vbarnew(i,:)] = KEP2RVmod(dkepdt(i,:));
end

% Plotting the changing orbit over time
figure(2)
rad = 6371;
[x,y,z] = sphere(1000);
set(gcf, 'Position', [300 300 1200 550])
mesh(x*rad,y*rad,z*rad); hold on;
plot3(rbarnew(:,1),rbarnew(:,2),rbarnew(:,3))




