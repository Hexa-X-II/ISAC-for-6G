function [FIM,d1, flag] = calculate_complete_fim_3d(params)
    % Initialize dimensions
    Nch = 7;  % Channel parameters per path [alpha_r, alpha_i, tau, phi, theta] for 3D
    flag = 0 ;
    n = 0:params.N - 1;
    % need to change to QAMCOM later
    [W, F] = generate_precoders_3d(params.Ntx_x, params.Ntx_y, params.Nrx_x, params.Nrx_y, 'DFT');

    k_indices = (0:floor(3300/params.K):floor(3300/params.K)*params.K-1).';
    exp_k = @(delay) exp(-1j*2*pi*params.delta_f * delay * k_indices);  % [K x 1] 

    reflected_paths = struct();
    reflected_paths_count = 0;

    % Wall reflections calculation
    for i = 1:length(params.walls)
        % Calculate virtual anchor (eq. 24) in 3D
        pVA = params.pTx - 2 * params.walls(i).n * ((params.pTx - params.walls(i).p)' * params.walls(i).n);

        % Calculate reflection point (eq. 25) in 3D
        numerator = (params.walls(i).p - pVA)' * params.walls(i).n;
        denominator = (params.pRx - pVA)' * params.walls(i).n;
        pr = pVA + (numerator/denominator) * (params.pRx - pVA);

        % Calculate vectors in local coordinate systems (3D)
        u_vec = params.RTx' * (pr - params.pTx);
        v_vec = params.RRx' * (pr - params.pRx);

        % Calculate angles for Tx (3D)
        u_norm_Tx = norm(u_vec);
        phi_Tx = atan2(u_vec(2), u_vec(1)) - pi/2;
        theta_Tx = asin(u_vec(3)/u_norm_Tx);

        % Calculate angles for Rx (3D)
        v_norm_Rx = norm(v_vec);
        phi_Rx = atan2(v_vec(2), v_vec(1)) - pi/2;
        theta_Rx = asin(v_vec(3)/v_norm_Rx);

        % Check both azimuth and elevation angles are within FoV
        if (abs(phi_Tx) < params.FoV_tx / 2) && (abs(theta_Tx) < params.FoV_tx_elevation / 2) && ...
           (abs(phi_Rx) < params.FoV_rx / 2) && (abs(theta_Rx) < params.FoV_rx_elevation / 2)

            reflected_paths_count = reflected_paths_count + 1;
            path_length = norm(u_vec) + norm(v_vec);
            tau = path_length / params.c;

            % 3D angle of incidence calculation (corrected)
            % Use the incident vector (pr - pTx) and the reflected vector (pRx - pr)
            angle_of_incidence = 0.5 * acos(((pr - params.pTx)' * (params.pRx - pr)) / ...
                                   (norm(pr - params.pTx) * norm(params.pRx - pr)));

            phase = -2 * pi * params.fc * path_length / params.c + params.phase_offset;
            gain = antenna_element_pattern(phi_Tx, theta_Tx, 65/180*pi, 8) * ...
                   antenna_element_pattern(phi_Rx, theta_Rx, 65/180*pi, 8);
            alpha = gain * params.wavelength * reflection_coeff(angle_of_incidence, params.epsilon_r) * ...
                    exp(1j * phase) / (4 * pi * path_length);

            % Store all 3D parameters
            reflected_paths(reflected_paths_count).alpha = alpha;
            reflected_paths(reflected_paths_count).tau = tau;
            reflected_paths(reflected_paths_count).phi_Rx = phi_Rx;
            reflected_paths(reflected_paths_count).phi_Tx = phi_Tx;
            reflected_paths(reflected_paths_count).theta_Rx = theta_Rx;
            reflected_paths(reflected_paths_count).theta_Tx = theta_Tx;
        end
    end

    L = 2 + length(reflected_paths);
    X = sqrt(params.Es) * ones(params.K, params.N);
    grad_all = zeros(L * Nch, 1, params.K, params.N);

    %% LOS path calculations (3D)
    u0 = params.RTx' * (params.pRx - params.pTx);
    v0 = params.RRx' * (params.pTx - params.pRx);

    % Calculate both azimuth and elevation angles
    phi_Tx0 = atan2(u0(2), u0(1)) - pi/2;
    theta_Tx0 = asin(u0(3) / norm(u0));
    phi_Rx0 = atan2(v0(2), v0(1)) - pi/2;
    theta_Rx0 = asin(v0(3) / norm(v0));

    % Use 3D steering vectors
    aTx0 = calculate_steeringvector_3d(phi_Tx0, theta_Tx0, params.QTx, params.wavelength);
    aRx0 = calculate_steeringvector_3d(phi_Rx0, theta_Rx0, params.QRx, params.wavelength);

    distance0 = norm(params.pRx - params.pTx);
    tau0 = distance0 / params.c;
    phase0 = -2 * pi * params.fc * distance0 / params.c + params.phase_offset;
    gain0 = antenna_element_pattern(phi_Tx0, theta_Tx0, 65/180*pi, 8) * ...
            antenna_element_pattern(phi_Rx0, theta_Rx0, 65/180*pi, 8);
    alpha0 = gain0 * params.wavelength * exp(1j * phase0) / (4 * pi * distance0);

    grad_all(1:Nch, 1, :, :) = calculate_gradient_matrix_3d(W, F, aTx0, aRx0, exp_k(tau0), ...
                                         X, alpha0, phi_Rx0, phi_Tx0, theta_Rx0, theta_Tx0, params);

   %function G = calculate_gradient_matrix_3d(W, F, aTx, aRx, exp_k, X, alpha, phi_rx, phi_tx, theta_rx, theta_tx, params)


    %% Target path calculations (3D)
    u1 = params.RTx' * (params.pT - params.pTx);
    v1 = params.RRx' * (params.pT - params.pRx);

    % Calculate both azimuth and elevation angles for target
    phi_Tx1 = atan2(u1(2), u1(1)) - pi/2;
    theta_Tx1 = asin(u1(3) / norm(u1));
    phi_Rx1 = atan2(v1(2), v1(1)) - pi/2;
    theta_Rx1 = asin(v1(3) / norm(v1));

    dTx_T = norm(params.pT - params.pTx);
    dT_Rx = norm(params.pRx - params.pT);
    d1=dTx_T + dT_Rx;
    tau1 = (dTx_T + dT_Rx) / params.c;

    % Use 3D steering vectors for target
    aTx1 = calculate_steeringvector_3d(phi_Tx1, theta_Tx1, params.QTx, params.wavelength);
    aRx1 = calculate_steeringvector_3d(phi_Rx1, theta_Rx1, params.QRx, params.wavelength);

    phase1 = -2 * pi * params.fc * (dTx_T + dT_Rx) / params.c + params.phase_offset + params.target_phase_offset;
    gain1 = antenna_element_pattern(phi_Tx1, theta_Tx1, 65/180*pi, 8) * ...
            antenna_element_pattern(phi_Rx1, theta_Rx1, 65/180*pi, 8);
    alpha1 = gain1 * params.wavelength * exp(1j * phase1) * params.sigma_RCS / ((4 * pi) ^ 1.5 * dTx_T * dT_Rx);

    % Check both azimuth and elevation angles for the target.
    % If the target is outside the FoV, terminate computation.
    if ~((abs(phi_Tx1) < params.FoV_tx / 2) && (abs(theta_Tx1) < params.FoV_tx_elevation / 2) && ...
        (abs(phi_Rx1) < params.FoV_rx / 2) && (abs(theta_Rx1) < params.FoV_rx_elevation / 2))
        
       disp('Target is outside FoV. Terminating computation.');
       flag = 1;
       % Here we return a default FIM.
       % The size is determined by the total number of paths, i.e., L = 2 + length(reflected_paths)
       % with Nch parameters per path.
       L = 2 + length(reflected_paths);
       FIM = zeros(L * Nch, L * Nch);
       return;
    end
    grad_all(Nch + 1:2 * Nch, 1, :, :) = calculate_gradient_matrix_3d(W, F, aTx1, aRx1, exp_k(tau1), ...
                                               X, alpha1, phi_Rx1, phi_Tx1, theta_Rx1, theta_Tx1, params);
    
    %% Wall reflections gradients
    count = 0;
    for ch_params = reflected_paths
        aTx = calculate_steeringvector_3d(ch_params.phi_Tx, ch_params.theta_Tx, params.QTx, params.wavelength);
        aRx = calculate_steeringvector_3d(ch_params.phi_Rx, ch_params.theta_Rx, params.QRx, params.wavelength);

        grad_all(2 * Nch + count * Nch + 1:2 * Nch + (count + 1) * Nch, 1, :, :) = ...
            calculate_gradient_matrix_3d(W, F, aTx, aRx, exp_k(ch_params.tau), X, ch_params.alpha, ...
                                         ch_params.phi_Rx, ch_params.phi_Tx, ch_params.theta_Rx, ch_params.theta_Tx, params);
        count = count + 1;
    end

    delay_Res = 1/(params.K*params.delta_f);
    angleTx_Res = 2/params.Ntx_x;
    angleRx_Res = 2/params.Ntx_x;

    delta_delay = abs(tau1-[tau0 reflected_paths.tau])<delay_Res;
    delta_angleTx = abs(sin(phi_Tx1)-[sin(phi_Tx0) sin([reflected_paths.phi_Tx])])<angleTx_Res;
    delta_angleRx = abs(sin(phi_Rx1)-[sin(phi_Rx0) sin([reflected_paths.phi_Rx])])<angleRx_Res;
    deltas = sum((delta_delay + delta_angleTx + delta_angleRx)==3,'all')>=1;


    % Compute FIM
    sigma2 = params.N0 * params.NNF;
    FIM = 2 / sigma2 * real(pagemtimes(grad_all, 'none', grad_all, 'ctranspose'));
    FIM = sum(FIM, [3, 4]);

    if deltas
        disp('Target cannot be resolved. Terminating computation.');
        flag = 2;
        % Here we return a default FIM.
        % The size is determined by the total number of paths, i.e., L = 2 + length(reflected_paths)
        % with Nch parameters per path.
        L = 2 + length(reflected_paths);
        % FIM = zeros(L * Nch, L * Nch);
        return;
     end
end

function gamma = reflection_coeff(theta_i, epsilon_r)
    % reflection_coeff Computes the Fresnel reflection coefficient for
    % perpendicular (s-polarized) waves.
    %
    %   The expression implemented is:
    %
    %       gamma = (sqrt(1/epsilon_r)*cos(theta_i) - sqrt(1 - (1/epsilon_r)*sin(theta_i).^2)) ...
    %               ./ (sqrt(1/epsilon_r)*cos(theta_i) + sqrt(1 - (1/epsilon_r)*sin(theta_i).^2));
    
        % Compute the terms in the numerator and denominator
        term1 = sqrt(1/epsilon_r) * cos(theta_i);
        term2 = sqrt(1 - (1/epsilon_r) * sin(theta_i).^2);
        
        % Compute the reflection coefficient
        gamma = (term1 - term2) ./ (term1 + term2);
    end
    