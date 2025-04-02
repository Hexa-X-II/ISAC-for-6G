function G = calculate_gradient_matrix_3d(W, F, aTx, aRx, exp_k, X, alpha, phi_rx, phi_tx, theta_rx, theta_tx, params)
% CALCULATE_GRADIENT_MATRIX_3D Computes the gradient matrix for a 3D 
% localization system.
%
%   G = calculate_gradient_matrix_3d(W, F, aTx, aRx, exp_k, X, alpha, 
%         phi_rx, phi_tx, theta_rx, theta_tx, params)

%
% Output:
%   G - A gradient matrix of size [7 x K x N] where 7 is the number of 
%       channel parameters per path.

    % Compute effective Tx and Rx contributions:
    ATx = repmat(aTx.' * F, params.K, 1);   % [K x N]
    ARx = repmat((W.' * aRx).', params.K, 1); % [K x N]
    B   = repmat(exp_k, 1, params.N);         % [K x N]
    
    % Combine contributions:
    C = ATx .* ARx .* B .* X;  % [K x N]
    
    % Initialize gradient matrix: 7 parameters per path.
    G = zeros(7, params.K, params.N);
    
    % Derivatives with respect to the real and imaginary parts of alpha:
    G(1, :, :) = C;         % d/d(Re{alpha})
    G(2, :, :) = 1j * C;      % d/d(Im{alpha})
    
    % Derivative with respect to delay (tau):
    % (0:K-1).' creates a [K x 1] column vector of subcarrier indices.
    G(3, :, :) = (-1j * 2 * pi * params.delta_f * (0:params.K - 1).') .* alpha .* C;
    
    % --- Derivatives with respect to Angle-of-Arrival (Rx side) ---
    % Compute the derivatives of the wave vector with respect to the angles,
    % passing the wavelength from params.
    [dkdphi_rx, dkdtheta_rx] = calculate_kvec_derivatives(phi_rx, theta_rx, params.wavelength);
    
    % Compute the derivative of the effective Rx steering vector.
    % Note: W.'*aRx gives an [N x 1] vector, so we transpose to get a row
    % vector and then replicate over K rows.
    dARxdphi   = repmat((W.' * (aRx .* (1j * (params.QRx.' * dkdphi_rx)))).', params.K, 1);
    dARxdtheta = repmat((W.' * (aRx .* (1j * (params.QRx.' * dkdtheta_rx)))).', params.K, 1);
    
    % Multiply element-wise to obtain the Rx derivatives.
    G(4, :, :) = alpha .* dARxdphi   .* ATx .* B .* X;  % derivative with respect to phi_rx
    G(6, :, :) = alpha .* dARxdtheta .* ATx .* B .* X;  % derivative with respect to theta_rx
    
    % --- Derivatives with respect to Angle-of-Departure (Tx side) ---
    [dkdphi_tx, dkdtheta_tx] = calculate_kvec_derivatives(phi_tx, theta_tx, params.wavelength);
    
    % For the Tx side, aTx.' * F yields a [1 x N] vector.
    dATxdphi   = repmat(( (aTx .* (1j * (params.QTx.' * dkdphi_tx))).' * F ), params.K, 1);
    dATxdtheta = repmat(( (aTx .* (1j * (params.QTx.' * dkdtheta_tx))).' * F ), params.K, 1);
    
    G(5, :, :) = alpha .* ARx .* dATxdphi   .* B .* X;  % derivative with respect to phi_tx
    G(7, :, :) = alpha .* ARx .* dATxdtheta .* B .* X;  % derivative with respect to theta_tx
end

function [dkdphi, dkdtheta] = calculate_kvec_derivatives(phi, theta, wavelength)
% CALCULATE_KVEC_DERIVATIVES Computes the derivative of the wave vector with
% respect to the azimuth (phi) and elevation (theta) angles.
%
%   [dkdphi, dkdtheta] = calculate_kvec_derivatives(phi, theta, wavelength)
%
% Outputs:
%   dkdphi    - Derivative of the wave vector with respect to phi (3x1).
%   dkdtheta  - Derivative of the wave vector with respect to theta (3x1).
%
% The derivatives are computed based on the standard definition of the 
% wave vector, \( k = \frac{2\pi}{\lambda}\hat{d} \).

    % Derivative with respect to the azimuth angle (phi):
    dkdphi = 2 * pi / wavelength .* [ -cos(phi)*cos(theta); 
                                        -sin(phi)*cos(theta); 
                                         0.0];
    % Derivative with respect to the elevation angle (theta):
    dkdtheta = 2 * pi / wavelength .* [ sin(phi)*sin(theta); 
                                        -cos(phi)*sin(theta); 
                                         cos(theta)];
end
