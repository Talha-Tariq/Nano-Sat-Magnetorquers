function  [x_dot] = ODEstest(t,x,output_flag,dummy_matix,t_span,waitbar_handle);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call constants
constants

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics

% First, extract states in a convenient form. 

eps = [x(1); x(2); x(3)];
%cross product operator
eps_x = crossop(eps);
eta = x(4);
omega = [x(5); x(6); x(7)];
omega_x = crossop(omega);

%{
r = [altitude*cos(2*pi/5400*t); 0*t; altitude*sin(2*pi/5400*t)];

persistent b 
b(:,end+1) = EarthMagField(r, t);  %MATLAB USES OWN TIMESTEPS :(
assignin('base','b_out',b);
%}

dot_eps = (0.5)*(eta*eye(3) + eps_x)*omega; 
dot_eta = (0.5)*(-transpose(eps)*omega);
torque_d = -c*omega - k*eps;
dot_omega = inv(inertia)*(-omega_x * inertia * omega + torque_d);
dot_x = [dot_eps; dot_eta; dot_omega];


%{
x1 = x(1);
x2 = x(2);

% Form dot_x = f(x,u) system.
dot_x1 = x2;
dot_x2 = 1/m*( - k1*x1 - k2*x1^3 - d*x2 + gamma1*cos(omega*t) );

dot_x = [dot_x1; dot_x2];

% Could also use this code
% dot_x = [x2; 1/m*( - k1*x1 - k2*x1^3 - d*x2 + gamma1*cos(omega*t) )];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output.
%}

if (output_flag == 0)
    % Output the rate of change of x.
    x_dot = dot_x;
    waitbar(t/t_span(length(t_span)),waitbar_handle);
    
elseif (output_flag == 1)   
    % This is for other kinds of post processing; not needed right now.
    x_dot = x;
    
end



