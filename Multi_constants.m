% Constants file.
% All constants should go in one m file. 

Re = 6371.2*1000; %radius of earth in m
altitude = 370*1000 + Re; %200,000 m from surface
orbit_period = 5400; %secs = 90 mins
%inertia = [0.1043, 0, 0; 0, 0.1020, 0; 0, 0, 0.0031]; %1kg mass 
inertia = [133/60000, 0, 0;0, 133/60000, 0;0, 0, 133/60000]; %1.33kg
dummy_matix = blkdiag([1 2; 3 4], 1);
I_min = min(inertia(inertia>0));
zeta = pi/2;
mean_motion = 2*pi/orbit_period; %rad/s
gain_op = 2*mean_motion*(1+sin(zeta))*I_min;
k = 0.00000001; %c^2/(4*norm(inertia)) %spring constant 
c = 0.00001; %damping constant 








