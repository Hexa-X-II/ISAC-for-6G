function a = calculate_steeringvector_3d(phi, theta, Q, lambda)
    % calculate_steeringvector_3d calculates the spatial steering vector for the 3D case.
    %
    %   Inputs:
    %     phi    - Azimuth angle (in radians)
    %     theta  - Elevation angle (in radians)
    %     Q      - Antenna position matrix. Expected to be 3 x N for a full 3D array.
    %              If Q is 2 x N (i.e., only x and y positions), a row of zeros is appended.
    %     lambda - Wavelength of the signal.
    %
    %   Output:
    %     a      - Spatial steering vector (N x 1).
    
        % If Q has only 2 rows, append a row of zeros for the z-coordinate.
        if size(Q, 1) == 2
            Q = [Q; zeros(1, size(Q, 2))];
        end
    
        % Calculate the 3D wave vector.
        % The formulation here is:
        %   k = (2*pi/lambda) * [ -sin(phi)*cos(theta);
        %                         cos(phi)*cos(theta);
        %                         sin(theta) ]
        k = (2 * pi / lambda) * [ -sin(phi)*cos(theta); 
                                   cos(phi)*cos(theta); 
                                   sin(theta) ];
        
        % Compute the steering vector.
        % Q.' is of size (N x 3) and k is (3 x 1), resulting in a (N x 1) vector.
        a = exp(1j * (Q.' * k));
    end
    