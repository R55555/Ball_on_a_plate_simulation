% Parameters

clear all; close all; clc;

% Ball Plate Parameters
m = 0.5;                % Mass of the ball (kg)
r = 0.05;               % Radius of the ball (m)
g = 9.81;               % Gravitational acceleration (m/s^2)
I = (2/5) * m * r^2;    % Moment of inertia of the ball (kg*m^2)

const_coeff = (m*g/(m + I/r^2)); % Constant in System Dynamics

% Trajectory parameters
r = 0.25;         % Radius of the trajectory (m)
omega = 0.50;     % Angular frequency (rad/s)


%%

%LQR

% Define linearised state-space model
A = [0, 1, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 1;
     0, 0, 0, 0];

B = [0,        0;
    -5/7 * g, 0;
     0,        0;
     0, -5/7 * g];

% Define LQR weighting matrices
Q = diag([10, 1, 10, 1]); % State weighting matrix
R = diag([1, 1]);         % Control input weighting matrix

% Calculate the LQR gain matrix
K_lqr = lqr(A, B, Q, R);


%%

% PID

% Time Step
dt = 0.01; 

% Setting the PID gains
K_pid = [0.2, 1., 0.25]; % KI, KP, KD


%%

% SMC

% Setting the SMC gains
% Sigma_K, Sigma_delta, S_n1, S_n2, 1 
K_smc = [2.5, 1.0, 3.85, 4.40, 1];


%%

% Define initial state and
% simulation duration
ti = 0; tf = 20;
X0 = [0.0; 0; 0; 0];
tspan = ti : dt : tf;

% Choosing the Trajectory
curr_traj = @(t) circularTrajectory(t, r, omega);
% curr_traj = @(t) fig8_Trajectory(t, r, omega);

% Choosing the respective controller
% xdot_ = @(t, x) xdot_lqr(t, x, const_coeff, K_lqr, curr_traj);
 xdot_ = @(t, x) xdot_pid(t, x, const_coeff, K_pid, dt, curr_traj);
%xdot_ = @(t, x) xdot_smc(t, x, const_coeff, K_smc, dt, curr_traj);

% Simulating the closed-loop system
U = zeros(length(tspan), 2);
X_actual = zeros(length(tspan), 4);
X_reference = zeros(length(tspan), 4);

for i = 1:length(tspan)
    [x_dot_curr, U(i,:), X_reference(i,:)] = xdot_(tspan(i), X_actual(i,:));
    X_actual(i+1,:) = X_actual(i,:) + x_dot_curr'*dt;
end

%%

% Inverse Kinematics
% D = 0.1;
% Rm = 0.03;
% rb = 0.75;
% r_pl = 0.1;
% rd  = 0.01;
% Tb = [0 0 0.07];

D = 0.23;
Rm = 0.08; 
rb = 0.17;
r_pl = 0.15;
%rd = 0.05;
rd = 0.0;
Tb= [0 0 0.22];

r= deg2rad([60 120 180 240 300 360]);
rt= deg2rad([90 90 240 240 300 300]);
rp= deg2rad([45 135 165 255 285 15]);
rt_ = deg2rad([180 0 300 120 60 240]);
%or  = [-1 +1 -1 +1 -1 +1]

B = zeros(6,3); 
B(:,1) = rb*cos(r)+ rd*cos(rt);
B(:,2) = rb*sin(r)+ rd*sin(rt);

P_p = zeros(6,3);
P_p(:,1) = r_pl*cos(rp);
P_p(:,2) = r_pl*sin(rp);

delta = zeros(length(tspan), 6);
modified_inv_kin = @(alpha, beta) inv_kinematics(alpha, beta, B, P_p, Tb, Rm, D, rt_);

for i = 1:length(tspan)
    delta(i,:) = modified_inv_kin(U(i,1), U(i,2))';
end


%%

% Plot results
figure;
hold on;
plot(X_actual(:,1),X_actual(:,3), 'r');
plot(X_reference(:,1), X_reference(:,3), 'b--');
hold off;

xlabel('X (m)'); ylabel('Y (m)');
legend('Actual', 'Reference');
axis equal; grid on;
title('Actual Vs Reference Trajectories');

% Plot results
figure;
plot(tspan, rad2deg(U));

xlabel('Time (s)'); ylabel('Control inputs (deg)');
legend('Alpha (deg)', 'Beta (deg)');
title('Plate Angles');
grid on;

% Plot results
figure;
error = X_reference - X_actual(1:end-1,:);
plot(tspan, error(:,1), tspan, error(:,3));

xlabel('Time (s)'); ylabel('Error');
legend('X (m)', 'Y (m)');
title('Error in Position');
grid on;

% Plot results
figure;
plot(tspan, delta)
xlabel('Time (s)'); ylabel('Join Angles (deg)');
title('Joint Angles Variation');
legend('Joint Angle 1', 'Joint Angle 2', 'Joint Angle 3', 'Joint Angle 4', 'Joint Angle 5', 'Joint Angle 6');
grid on;
%% Calculating performance metrics for X and Y positions

% Calculate error in position (only X and Y positions are relevant)
error_x = error(:,1); % X-error over time
error_y = error(:,3); % Y-error over time

% Settling time (within 2% tolerance)
tolerance = 0.05;
% Define a minimum tolerance threshold
%min_tolerance_threshold = 0.01;  % for example, 1 cm

% Calculate tolerance as the max of (percentage-based tolerance OR minimum threshold)
settling_idx_x = find(abs(error_x) > tolerance*max(X_reference(:,1)), 1, 'last');
settling_idx_y = find(abs(error_y) > tolerance*max(X_reference(:,3)), 1, 'last');
settling_time_x = tspan(settling_idx_x);
settling_time_y = tspan(settling_idx_y);
settling_time = max(settling_time_x, settling_time_y)

 % Rise time (from 10% to 90% of the final value)
 rise_time_x = find(abs(X_actual(1:end-1,1)) >= 0.9 * abs(X_reference(:,1)), 1, 'first') * dt;
 rise_time_y = find(abs(X_actual(1:end-1,3)) >= 0.9 * abs(X_reference(:,3)), 1, 'first') * dt;
 rise_time = max(rise_time_x, rise_time_y)
% 
% % Find the time when actual position first reaches 10% and 90% of reference position
% rise_time_x_start = find(abs(X_actual(:,1)) >= 0.1 * abs(X_reference(:,1)), 1, 'first') * dt;
% rise_time_x_end = find(abs(X_actual(:,1)) >= 0.9 * abs(X_reference(:,1)), 1, 'first') * dt;
% 
% rise_time_y_start = find(abs(X_actual(:,3)) >= 0.1 * abs(X_reference(:,3)), 1, 'first') * dt;
% rise_time_y_end = find(abs(X_actual(:,3)) >= 0.9 * abs(X_reference(:,3)), 1, 'first') * dt;
% 
% % Calculate rise time as the time between reaching 10% and 90% for both axes
% rise_time_x = rise_time_x_end - rise_time_x_start;
% rise_time_y = rise_time_y_end - rise_time_y_start;

% Overall rise time
%rise_time = max(rise_time_x, rise_time_y);


% Steady-state error calculation
steady_state_error_x = max(X_reference(settling_idx_x:end,1) - X_actual(settling_idx_x:end-1,1));
steady_state_error_y = max(X_reference(settling_idx_x:end,3) - X_actual(settling_idx_x:end-1,3));
steady_state_error = max(steady_state_error_x, steady_state_error_y)
