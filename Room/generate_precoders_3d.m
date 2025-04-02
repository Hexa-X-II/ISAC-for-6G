function [W, F] = generate_precoders_3d(Ntx_x, Ntx_z, Nrx_x, Nrx_z, type)
% [W, F] = generate_precoders_3d(Ntx_x, Ntx_z, Nrx_x, Nrx_z)
%   Generate DFT precoders for the 3D case. Arrays are Ntx_z x Ntx_x (Tx)
%   and Nrx_z x Nrx_z (Rx) in size.
%
%   Inputs:
%   -   Ntx_x:  Number of Tx elements along x.
%   -   Ntx_z:  Number of Tx elements along z.
%   -   Nrx_x:  Number of Rx elements along x.
%   -   Nrx_z:  Number of Rx elements along z.
%
%   Outputs:
%   -   W:      Rx precoder sequence, (Nrx_x * Nrx_z) x N
%   -   F:      Tx precider sequence, (Ntx_x * Ntx_z) x N

Nrx = Nrx_x * Nrx_z;
Ntx = Ntx_x * Ntx_z;

xs_rx = 0:Nrx_x - 1;
zs_rx = (0:Nrx_z - 1).';

xs_tx = 0:Ntx_x - 1;
zs_tx = (0:Ntx_z - 1).';

switch type
    case "DFT"
        N = Ntx_x * Ntx_z * Nrx_x * Nrx_z;
        
        W = zeros(Nrx, N);
        F = zeros(Ntx, N);
        for n = 1:N
            
            % See doc
            nurx = mod(mod(n, Nrx), Nrx_x) - (Nrx_x - 1) / 2;
            gammarx = floor(mod(n, Nrx) / Nrx_x) - (Nrx_z - 1) / 2;
        
            nutx = mod(floor(n / Nrx), Ntx_x) - (Ntx_x - 1) / 2;
            gammatx = floor(floor(n / Nrx) / Ntx_x) - (Ntx_z - 1) / 2;
        
            w = 1 / sqrt(Nrx) .* exp(-1j .* 2 .* pi .* (xs_rx .* nurx / Nrx_x + zs_rx .* gammarx / Nrx_z));
            f = 1 / sqrt(Ntx) .* exp(-1j .* 2 .* pi .* (xs_tx .* nutx / Ntx_x + zs_tx .* gammatx / Ntx_z));
        
            W(:, n) = w(:);
            F(:, n) = f(:);
        end
    case "QAMCOM"

        span_tx = 90 / 180 * pi;
        span_rx = 90 / 180 * pi;
        num_beams_tx = 10; % Number of beams, Tx
        num_beams_rx = 10; % Number of beams, Rx
        N = num_beams_tx  * num_beams_rx;

        % Sweeping determined by angular span and number of beams
        n = 0:N - 1;
        nPrime = floor(n / num_beams_rx);
        nBis = mod(n, num_beams_rx);

        anglesTx = span_tx / num_beams_tx * (nPrime - (num_beams_tx - 1) / 2);
        anglesRx = span_rx / num_beams_rx * (nBis - (num_beams_rx - 1) / 2);

        QTx = generate_Q_matrix(Ntx_x, Ntx_z, 1.0); % Using bogus wavelength, doesn't matter as long as I'm consistent with the kvec
        QRx = generate_Q_matrix(Nrx_x, Nrx_z, 1.0);

        W = zeros(Nrx, N);
        F = zeros(Ntx, N);

        for i = 1:N
            kvec_tx = compute_kvec(anglesTx(i), 0.0, 1.0); % steering only in azimuth (phi), bogus wavelength
            F(:, i) = 1 / sqrt(Ntx) .* exp(-1j .* QTx.' * kvec_tx);

            kvec_rx = compute_kvec(anglesRx(i), 0.0, 1.0); 
            W(:, i) = 1 / sqrt(Nrx) .* exp(-1j .* QRx.' * kvec_rx);
        end
end
end

function kvec = compute_kvec(phi, theta, lambda)
    kvec = [-sin(phi) .* cos(theta); ...
            cos(phi) .* cos(theta); ...
            sin(theta)] * 2 * pi / lambda; 
end