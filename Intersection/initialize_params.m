%% Helper Functions
function params = initialize_params(intersection_size)
    params.c = 299792458;
    params.fc = 10e9;
    
    switch params.fc
        case 10e9
            params.delta_f = 120e3;
            params.Ntx_x = 2;
            params.Ntx_y = 2;
            params.Nrx_x = 2;
            params.Nrx_y = 2;
        case 60e9
            params.delta_f = 960e3;
            params.Ntx_x = 4;
            params.Ntx_y = 4;
            params.Nrx_x = 4;
            params.Nrx_y = 4;
        case 100e9
            params.delta_f = 1920e3;
            params.Ntx_x = 8;
            params.Ntx_y = 4;
            params.Nrx_x = 4;
            params.Nrx_y = 4;
    end
    params.K = 792;
    params.wavelength = params.c / params.fc;
    params.BW = params.delta_f * params.K;
    
    params.Es = 60;    
    params.Es = 10^(params.Es/10)*1e-3;
    params.Es = params.Es / params.BW;
    
    params.kB = 1.380649e-23;
    params.T0 = 290;
    params.N0 = params.kB * params.T0;
    params.NNF = 10^(7/10);
    
    params.Ntx = params.Ntx_x * params.Ntx_y;
    params.Nrx = params.Nrx_x * params.Nrx_y;
    
    params.QTx = generate_Q_matrix(params.Ntx_x, params.Ntx_y, params.wavelength);
    params.QRx = generate_Q_matrix(params.Nrx_x, params.Nrx_y, params.wavelength);
    params.N = params.Ntx * params.Nrx;
    
    % Place Tx on the top street, 1 m above ground.
    params.pTx = [25; 49; 1];
    
    % Place Rx on top of the bottom-right building.
    alpha = 0.3;
    ix_end   = intersection_size(1) * (1 - alpha);  
    iy_start = intersection_size(2) * alpha;        
    params.pRx = [ix_end; iy_start; intersection_size(3) + 1];  % [35; 15; 6]
    
    params.pT  = [2; 2; 1];
    % theta is the pitch (up and down along the x axis)
    % phi is the yaw/azimith left and right along the z axis
    % psi is the roll along the y axis
    % all counter clockwise
    params.TxTheta = 0;
    params.TxPhi   = pi;
    params.TxPsi   = 0;
    % looking down for 30 degree
    params.RxTheta = -pi/6;
    params.RxPhi   = pi/4;
    params.RxPsi   = 0;
    params.RTx = euler_rotation_matrix(params.TxTheta, params.TxPhi, params.TxPsi);
    params.RRx = euler_rotation_matrix(params.RxTheta, params.RxPhi, params.RxPsi);
    
    params.FoV_tx = pi - 0.0001;
    params.FoV_rx = pi - 0.0001;
    params.FoV_tx_elevation = pi - 0.0001;
    params.FoV_rx_elevation = pi - 0.0001;
    
    alpha = 0.3;
    ix_start = intersection_size(1)*alpha;
    ix_end   = intersection_size(1)*(1 - alpha);
    iy_start = intersection_size(2)*alpha;
    iy_end   = intersection_size(2)*(1 - alpha);
    
    z_mid = intersection_size(3)/2;
    
    % this is the ground
    % check what is the ground, plot it.
    params.walls(1).p = [intersection_size(1)/2; intersection_size(2)/2; 0];
    params.walls(1).n = [0; 0; 1];
    params.walls(1).start = [0; 0; 0];
    params.walls(1).end   = [intersection_size(1); intersection_size(2); 0];
    
    % bottom left verticle (BLV)
    params.walls(2).p = [ix_start; iy_start/2; z_mid];
    params.walls(2).n = [1; 0; 0];
    params.walls(2).start = [ix_start; 0; z_mid];
    params.walls(2).end   = [ix_start; iy_start; z_mid];
    %BLH
    params.walls(3).p = [ix_start/2; iy_start; z_mid];
    params.walls(3).n = [0; 1; 0];
    params.walls(3).start = [0; iy_start; z_mid];
    params.walls(3).end   = [ix_start; iy_start; z_mid];
    % BRV
    params.walls(4).p = [ix_end; iy_start/2; z_mid];
    params.walls(4).n = [-1; 0; 0];
    params.walls(4).start = [ix_end; 0; z_mid];
    params.walls(4).end   = [ix_end; iy_start; z_mid];
    % BRH
    params.walls(5).p = [(ix_end + intersection_size(1))/2; iy_start; z_mid];
    params.walls(5).n = [0; 1; 0];
    params.walls(5).start = [ix_end; iy_start; z_mid];
    params.walls(5).end   = [intersection_size(1); iy_start; z_mid];
    % TLV
    params.walls(6).p = [ix_start; (iy_end + intersection_size(2))/2; z_mid];
    params.walls(6).n = [1; 0; 0];
    params.walls(6).start = [ix_start; iy_end; z_mid];
    params.walls(6).end   = [ix_start; intersection_size(2); z_mid];
    % TLH
    params.walls(7).p = [ix_start/2; iy_end; z_mid];
    params.walls(7).n = [0; -1; 0];
    params.walls(7).start = [0; iy_end; z_mid];
    params.walls(7).end   = [ix_start; iy_end; z_mid];
    % TRV
    params.walls(8).p = [ix_end; (iy_end + intersection_size(2))/2; z_mid];
    params.walls(8).n = [-1; 0; 0];
    params.walls(8).start = [ix_end; iy_end; z_mid];
    params.walls(8).end   = [ix_end; intersection_size(2); z_mid];
    % TRH
    params.walls(9).p = [(ix_end+intersection_size(1))/2; iy_end; z_mid];
    params.walls(9).n = [0; -1; 0];
    params.walls(9).start = [ix_end; iy_end; z_mid];
    params.walls(9).end   = [intersection_size(1); iy_end; z_mid];
    
    params.sigma_RCS = 1.0;
    params.epsilon_r = 5.0;
    params.phase_offset = pi/3;
    params.target_phase_offset = pi/8;
end

function R = euler_rotation_matrix(theta, phi, psi)
    Rx = [1, 0, 0;
          0, cos(theta), -sin(theta);
          0, sin(theta), cos(theta)];
    Rz = [cos(phi), -sin(phi), 0;
          sin(phi), cos(phi), 0;
          0, 0, 1];
    Ry = [cos(psi), 0, sin(psi);
          0, 1, 0;
          -sin(psi), 0, cos(psi)];
    R = Ry * Rx * Rz;
end