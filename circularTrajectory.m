function xref = circularTrajectory(t, r, omega)
    % Function to calculate the reference trajectory for a circle
    % Inputs:
    %   t - Time (s)
    %   r - Radius of the circle (m)
    %   omega - Angular frequency (rad/s)
    % Outputs:
    %   xref(1) - Reference x position at time t
    %   xref(2) - Reference x veloctiy at time t
    %   xref(3) - Reference y position at time t
    %   xref(4) - Reference y veloctiy at time t
    
    xref = zeros(1, 4);
    xref(1) = r * cos(omega * t);
    xref(3) = r * sin(omega * t);
    xref(2) = -r * omega * sin(omega * t);
    xref(4) = r * omega * cos(omega * t);
    
end
