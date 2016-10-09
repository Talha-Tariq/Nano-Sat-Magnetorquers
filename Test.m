
constants;

torque = (-c*x_out(:,5:7) - k*x_out(:,1:3)).';

%C = ((2*x_out(:,4).^2 -1).'*eye(t_div)).' + 2*x_out(:,1:3)*transpose(x_out(:,1:3)) - 2*x_out(:,4)*crossop(x_out(:,1:3));

%first iteration
C1 = (2*x_out(1,4)^2 - 1)*eye(3) + 2*x_out(1,1:3).'*(x_out(1,1:3)) - 2*x_out(1,4)*crossop(x_out(1,1:3));

%second iteration
C2 = (2*x_out(2,4)^2 - 1)*eye(3) + 2*x_out(2,1:3).'*(x_out(2,1:3)) - 2*x_out(2,4)*crossop(x_out(2,1:3));

