
clear all
clc
close all
format long

% profile on; % Try uncommenting this ``profile off;'' code at the bottom.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants 
constants % All constants in one file. 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions.
IC = [0 0 0 1 0.1 0.1 0.1].'; 

%IC = [0 0 0 1 0 -2*pi/orbit_period 0].'; %erase torque from dynamics eqn

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time.
t0 = 0; % s
t_max = 100000; % s  %60
t_div = 120001;    %1001
t_span = linspace(t0,t_max,t_div); % Total simulation time.
% t_span = [0 t_max];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation options.
% options = odeset('AbsTol',1e-10,'RelTol',1e-10); % This changes the integration tolerence.
options = [];

tic
waitbar_handle = waitbar(0,'Numerical Integration Progress');
output_flag = 0;

% Use any solver; I use these three the most. 
[t,x_out] = ode45(@ODEs,t_span,IC,options,output_flag,dummy_matix,t_span,waitbar_handle);
% [t,x_out] = ode15s(@ODEs,t_span,IC,options,output_flag,dummy_matix,t_span,waitbar_handle);
% [t,x_out] = ode113(@ODEs,t_span,IC,options,output_flag,dummy_matix,t_span,waitbar_handle);

%[t,x_out] = ode45(@Test,t_span,IC,options,output_flag,dummy_matix,t_span,waitbar_handle);

close(waitbar_handle);
time_stamp = toc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing
%post_processing

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot omega and epsilon data vs. time
plot_script_v1

% Check if f() = mag(epsilon)^2 + mag(eta)^2 - 1 = 0 to see if code is fine
accuracy_check

% Check if energy is constant
energy_check

MagMatrix

%Test

% profile off;
% profview;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save all the data. (You never know when you'll need it again.)
save sim_data_v1