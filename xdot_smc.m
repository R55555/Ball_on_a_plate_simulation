function [xdot, u, x_desired] = xdot_smc(curr_time, x_current, const_coeff, gain_K, dt, traj)
    global smc_error;

    % Compute the reference trajectory and Control Input
    x_desired = traj(curr_time);

    if(curr_time>0)
        smc_error(2:end,1) = x_current(1:2) - x_desired(1:2);
        smc_error(2:end,2) = x_current(3:end) - x_desired(3:end);
        smc_error(1,:) = smc_error(1,:) + smc_error(2,:)*dt;
    else
        smc_error = zeros(3, 2);
    end
    
    S = gain_K(3:end) * smc_error;
    u = (gain_K(3:end-1) * smc_error(2:end,:) + (gain_K(1) * S./(gain_K(2) + abs(S))))./const_coeff;
    
    % Compute derivatives
    xdot = zeros(4, 1);
    xdot(1) = x_current(2);                 % x dot
    xdot(2) = -const_coeff * sin(u(1));     % x double dot
    xdot(3) = x_current(4);                 % y dot
    xdot(4) = -const_coeff * sin(u(2));     % y double dot
end