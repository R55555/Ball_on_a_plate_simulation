function delta = inv_kinematics(alpha, beta, B, P_p, Tb, Rm, D, rt_)
    % Rotation matrices for alpha and beta
    Rx = [1 0 0; 0 cos(beta) sin(beta); 0 -sin(beta) cos(beta)];
    Ry = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
    Rb = Ry * Rx;
    
    % Transform platform points
    P = P_p * Rb' + Tb;
    
    % Ensure P and B have compatible dimensions
    % if size(P, 2) ~= 6 || size(B, 2) ~= 6
    %     error('P and B must each have 6 columns.');
    % end
    
    % Calculate vector differences
    L = P - B;
    delta = zeros(1, 6);

    % Calculate delta for each column of L
    for i = 1:6
        Li_norm = vecnorm(L(i, :))
        % Clamp argument of acos within [-1, 1] to avoid complex numbers
        asin_arg = (Rm^2 + Li_norm^2 - D^2) / (2 * Rm * Li_norm);
        %acos_arg = min(max(acos_arg, -1), 1);
        delta(i) = asin(asin_arg);
    end
    
    % Convert delta to degrees
    delta = rad2deg(delta);
end
