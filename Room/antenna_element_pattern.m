function g = antenna_element_pattern(phi, theta, hpbw, gmax)
% g = antenna_element_pattern(phi, theta, hpbw, gmax)
%   Cosine-shaped antenna element pattern function for 3D.
%
%   Inputs:
%   -   phi:    Azimuth angle sample point.   
%   -   theta   Elevation angle sample point.
%   -   hpbw:   Half-power beamwidth of antenna element pattern.
%   -   gmax:   Maximum gain (directivity) in dBi.
%
%   Outputs:
%   -   g:      Antenna amplitude pattern value at (phi, theta), 
%               dimension sqrt(Power).

exponent = - 0.5 * log(2) / log(cos(hpbw / 2));
smallest_level = 1e-3;

g = 10 .^ (gmax / 20) .* ... 
    abs(cos(phi)) .^ exponent .* ...
    abs(cos(theta)) .^ exponent + smallest_level;
end