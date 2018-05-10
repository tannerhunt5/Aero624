function dkepdt = ODEJ2Gauss(t, Kep)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global mu J2 Re

h     = Kep(1);
e     = Kep(2);
Omega = Kep(3);
inc   = Kep(4);
argp  = Kep(5);
TA    = Kep(6);

r = h^2/(mu*(1 + e*cos(TA)));
u = argp + TA;

% Solving for changing elements
const = (3/2)*((J2*mu*Re^2)/(r^3));
dhdt      = -3/2*J2*mu*Re^2/r^3*sin(inc)^2*sin(2*u);
dedt      = 3/2*J2*mu*Re^2/h/r^3*(h^2/mu/r ...
            *(sin(u)*sin(inc)^2*(3*sin(TA)*sin(u) - 2*cos(TA)*cos(u)) - sin(TA)) ...
            -sin(inc)^2*sin(2*u)*(e + cos(TA)));
dOmegadt  = -3*((J2*mu*Re^2)/(h*r^3))*(sin(u)^2*cos(inc));
dincdt    = -3/4*J2*mu*Re^2/h/r^3*sin(2*u)*sin(2*inc);
dargpdt   = const/e/h*(h/mu/r*cos(TA)*(1-3*sin(inc)^2*sin(u)^2) - (2 + e*cos(TA))*sin(2*u)*sin(inc)^2*sin(TA) + 2*e*cos(inc)^2*sin(u)^2);
dTAdt     = h/r^2 + const/e/h * (h^2/(mu*r)*cos(TA)*(3*sin(inc)^2*sin(u)^2 - 1) + (2 + e*cos(TA))*sin(2*u)*sin(inc)^2*sin(TA));

dkepdt = [dhdt dedt dOmegadt dincdt dargpdt dTAdt]';

end

