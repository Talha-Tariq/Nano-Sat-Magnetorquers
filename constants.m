% Constants file.
% All constants should go in one m file. 

%{
m = 1.5; % kg, mass.
k1 = 5; % N/m, spring constant.
k2 = 0.5; % N/m, second spring constant.
d = 0.35; % N/(m/2), damping constant.
gamma1 = 0.1*0; % N, amplitute of forcing function; note ``gamma'' is a function, which is why I use ``gamma1'' as the variable. Try ``help gamma'' in the command window.
omega = 10;; % rad/s, frequency of forcing function.
%}

Re = 6371.2*1000; %radius of earth in m
altitude = 200*1000 + Re; %200,000 m from surface
inertia = [1/600, 0, 0;0, 1/600, 0;0, 0, 1/600]; %[5, 0, 0;0, 10, 0;0, 0, 15]
dummy_matix = blkdiag([1 2; 3 4], 1);
k = 1; %c^2/(4*norm(inertia)) %spring constant %0.75
c = 0.01; %damping constant %1.5