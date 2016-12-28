function  [x_dot] = PD_ODEs(t,x,output_flag,dummy_matix,t_span,waitbar_handle);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call constants
PD_constants

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics

% First, extract states in a convenient form. 

eps = [x(1); x(2); x(3)];
%cross product operator
eps_x = crossop(eps);
eta = x(4);
omega = [x(5); x(6); x(7)];
omega_x = crossop(omega);

dot_eps = (0.5)*(eta*eye(3) + eps_x)*omega; 
dot_eta = (0.5)*(-transpose(eps)*omega);
torque_d = -c*omega - k*eps;


%time taken to reach quaternion tolerance
persistent t_quaternian
if (t > 5000) && (eps(1) <= 0.01) && (eps(2) <= 0.01) && (eps(3) <= 0.01) && (norm(t_quaternian) == 0)
    t_quaternian = t;
    assignin('base','t_quaternian',t_quaternian);
end

%PD control
r = [altitude*cos((2*pi/orbit_period)*t), 0*t, altitude*sin((2*pi/orbit_period)*t)];

b_eci = EarthMagField(r.', t);

rotation_mat = (2*(eta^2) - 1)*eye(3) + 2*eps*(eps.') - 2*eta*crossop(eps);

b_body = rotation_mat * b_eci;

torque_t = torque_d - (dot(torque_d,b_body)*b_body)/norm(b_body)^2;

dot_omega = inv(inertia)*(-omega_x * inertia * omega + torque_t);
dot_x = [dot_eps; dot_eta; dot_omega];

if (output_flag == 0)
    % Output the rate of change of x.
    x_dot = dot_x;
    waitbar(t/t_span(length(t_span)),waitbar_handle);
    
elseif (output_flag == 1)   
    % This is for other kinds of post processing; not needed right now.
    x_dot = x;
    
end



