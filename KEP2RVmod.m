function [rbar, vbar] = KEP2RVmod(Kep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
global mu

h     = Kep(1);
e     = Kep(2);
Omega = Kep(3);
inc   = Kep(4);
argp  = Kep(5);
TA    = Kep(6);

rp =(h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

R3Omega = [ cos(Omega) sin(Omega) 0;
        -sin(Omega) cos(Omega) 0;
        0 0 1];
    
R1i = [ 1 0 0;
        0 cos(inc) sin(inc);
        0 -sin(inc) cos(inc)];
    
R3argp = [ cos(argp) sin(argp) 0;
        -sin(argp) cos(argp) 0;
        0 0 1];

RotMat = (R3argp*R1i*R3Omega)';

r = RotMat*rp;
v = RotMat*vp;

rbar = r';
vbar = v';

end

