function xref= fig8_Trajectory(t, a, omega)
    % Function to calculate the reference trajectory for a circle
    % Inputs:
    %   t - Time (s)
    %   a - c*sqrt(2) 2c is foci (m)
    %   omega - Angular frequency (rad/s)
    % Outputs:
    %   xref(1) - Reference x position at time t
    %   xref(2) - Reference x veloctiy at time t
    %   xref(3) - Reference y position at time t
    %   xref(4) - Reference y veloctiy at time t

    xref = zeros(1, 4);
    xref(1) = a * sin(omega*t);
    xref(2) = a*omega*cos(omega * t);
    xref(3) = a/2 * sin(2*omega*t);
    xref(4) = a*omega*cos(2*omega * t);    
end
