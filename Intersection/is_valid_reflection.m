function valid = is_valid_reflection(Tx, Rx, wall)
    % IS_VALID_REFLECTION determines whether a specular reflection from Tx to Rx
    % off the specified wall is valid.
    %
    %   valid = is_valid_reflection(Tx, Rx, wall) returns true if the reflection
    %   is valid and false otherwise.
    %
    %   Inputs:
    %       Tx   - Transmitter position (vector)
    %       Rx   - Receiver position (vector)
    %       wall - Structure with the following fields:
    %                p     - A point on the wall plane
    %                n     - The unit normal vector of the wall
    %                start - Start point of the wall segment
    %                end   - End point of the wall segment
    %
    %   The function works by:
    %       1. Reflecting Tx across the wall plane.
    %       2. Determining the intersection point (R_point) between the line from
    %          the reflected point to Rx and the wall plane.
    %       3. Checking that this intersection point lies on the wall segment.
    
    % Reflect Tx across the plane defined by wall.p and wall.n
    Tx_ref = Tx - 2 * dot(Tx - wall.p, wall.n) * wall.n;
    
    % Calculate the direction from the reflected Tx to Rx
    d = Rx - Tx_ref;
    
    % Compute the dot product of d with the wall normal (denom for the ray-plane intersection)
    d_dot_n = dot(d, wall.n);
    
    % If d_dot_n is very close to zero, the ray is nearly parallel to the wall plane
    if abs(d_dot_n) < 1e-6
        valid = false;
        return;
    end
    
    % Compute the parameter t for the intersection point along the ray
    t = -dot(Tx_ref - wall.p, wall.n) / d_dot_n;
    
    % The intersection must occur between Tx_ref and Rx
    if t < 0 || t > 1
        valid = false;
        return;
    end
    
    % Determine the reflection point on the wall plane
    R_point = Tx_ref + t * d;
    
    % Check if R_point lies within the wall segment
    d_wall = wall.end - wall.start;
    proj = dot(R_point - wall.start, d_wall) / (norm(d_wall)^2);
    valid = (proj >= 0) && (proj <= 1);
end
