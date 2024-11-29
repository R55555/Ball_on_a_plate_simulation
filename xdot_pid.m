function [xdot, u, x_desired] = xdot_pid(curr_time, x_current, const_coeff, gain_K, dt, traj)
    global pid_error;
    
    % Compute the reference trajectory and Control Input
    x_desired = traj(curr_time);

    if(curr_time>0)
        pid_error(2:end,1) = x_desired(1:2) - x_current(1:2);
        pid_error(2:end,2) = x_desired(3:end) - x_current(3:end);
        pid_error(1,:) = pid_error(1,:) + pid_error(2,:)*dt;
    else
        pid_error = zeros(3, 2);
    end

    u = -gain_K * pid_error;

    % Compute derivatives
    xdot = zeros(4, 1);
    xdot(1) = x_current(2);                 % x dot
    xdot(2) = -const_coeff * sin(u(1));     % x double dot
    xdot(3) = x_current(4);                 % y dot
    xdot(4) = -const_coeff * sin(u(2));     % y double dot
end
