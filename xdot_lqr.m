function [xdot, u, x_desired] = xdot_lqr(curr_time, x_current, const_coeff, gain_K, traj)
    % Compute the reference trajectory and Control Input
    x_desired = traj(curr_time);
    u = -gain_K * (x_current - x_desired)';
    
    % Compute derivatives
    xdot = zeros(4, 1);
    xdot(1) = x_current(2);                 % x dot
    xdot(2) = -const_coeff * sin(u(1));     % x double dot
    xdot(3) = x_current(4);                 % y dot
    xdot(4) = -const_coeff * sin(u(2));     % y double dot
    u = u';
end
