function [TPEB, CRB_xT, CRB_yT, CRB_zT, flag] = calculate_position_error_bound_3d(FIM, params, flag)
    % Calculate Target Position Error Bound following Section 3 (3D version)
    % Input: FIM    - Full Fisher Information Matrix
    %        params - Parameters structure
    %
    % Step 1: Treat reflected paths as nuisance parameters.
    % From equations (63-64):
    % E(η') = B(η',η') - B(η',ξ) B^(-1)(ξ,ξ) B^T(η',ξ)
    %
    % Nch: Number of channel parameters per path.
    % For 3D, we have 7 parameters per path:
    %   [Re{α}, Im{α}, τ, φ_Rx, φ_Tx, θ_Rx, θ_Tx]
    Nch = 7;
    
    % Extract blocks from FIM following equation (64)
    % First 14 parameters are for LOS and target paths (2 paths * 7 parameters)
    B_eta_eta = FIM(1:2*Nch, 1:2*Nch);
    B_eta_xi  = FIM(1:2*Nch, 2*Nch+1:end);
    B_xi_xi   = FIM(2*Nch+1:end, 2*Nch+1:end);
    
    % Calculate equivalent FIM for LOS and target paths using equation (63)
    E_eta = B_eta_eta - B_eta_xi * inv(B_xi_xi) * B_eta_xi.';
    
    % Step 2: Transform parameters following equations (65-66)
    % Transform physical parameter vector:
    %   η'' = [α₀ (2); B (1); φ₀ (2); θ₀ (2); α₁ (2); p_T (3)]
    % (i.e., a 12×1 vector) to the full channel parameter vector.
    J = calculate_transformation_jacobian_3d(params);
    E_transformed = J.' * E_eta * J;
    
    % Step 3: Calculate TPEB using equation (73)
    CRB_matrix = inv(E_transformed);
    position_CRB = CRB_matrix(10:12, 10:12);
    
    % Calculate individual CRBs for x, y, and z coordinates.
    CRB_xT = position_CRB(1,1);
    CRB_yT = position_CRB(2,2);
    CRB_zT = position_CRB(3,3);
    
    
    % Calculate TPEB using equation (73)
    TPEB = sqrt(trace(position_CRB));
    if isnan(TPEB)
        flag = 1;
    end
end

function J = calculate_transformation_jacobian_3d(params)
% calculate_transformation_jacobian_3d computes the Jacobian of the mapping
% from the “physical” (or transformed) parameter vector
%
%      η'' = [α₀ (2); B (1); φ₀ (2); θ₀ (2); α₁ (2); p_T (3)]
%
% (i.e., a 12×1 vector) to the full channel parameter vector
%
%      η = [η₀; η₁] with η₀ ∈ ℝ⁷ (LOS) and η₁ ∈ ℝ⁷ (target),
%      so that η ∈ ℝ¹⁴.
%
% By design we write the overall Jacobian as
%
%      J_full = [ I_{7×7}         0_{7×5} ;
%                 J_target (7×5)             ],
%
% where J_target is defined as
%
%      J_target = [ I_{2×2}    0_{2×3} ;
%                   0_{5×2}  (∂g₁/∂p_T) (5×3) ],
%
% with
%      g₁ = [τ₁; φ_Rx,₁; φ_Tx,₁; θ_Rx,₁; θ_Tx,₁] ∈ ℝ⁵,
% i.e. the target’s delay and angular parameters expressed as functions
% of the target position p_T.
%
% Input:
%   params - structure that must include at least:
%            .pT   : target position [x_T; y_T; z_T]
%            .pTx  : transmitter position
%            .pRx  : receiver position
%            .RTx  : Tx rotation matrix (3×3)
%            .RRx  : Rx rotation matrix (3×3)
%            .c    : speed of light
%
%            (Other fields such as the estimated elevation angles may be
%             computed here from the geometry.)
%
% Output:
%   J - 14×12 transformation Jacobian

    % Initialize full Jacobian (14×12)
    J = zeros(14,12);
    
    % Upper block (rows 1:7): LOS channel parameters are assumed to be
    % directly taken from the first 7 entries of the transformed vector.
    % LOS physical parameters: [α₀ (2); B (1); φ₀ (2); θ₀ (2)]
    J(1:7, 1:7) = eye(7);
    % (Columns 8–12 remain zero for the LOS block.)
    
    % Lower block (rows 8:14) corresponds to the target path.
    % The target physical parameters are composed of:
    %   α₁ (2) and p_T (3)  → total 5 parameters.
    % We assume that the target channel parameters are given by:
    %   [α₁; g₁] ∈ ℝ⁷, where g₁ = [τ₁; φ_Rx,₁; φ_Tx,₁; θ_Rx,₁; θ_Tx,₁].
    %
    % For the first 2 entries (α₁), we assume a direct mapping.
    J(8:9, 8:9) = eye(2);
    % The remaining 5 entries (rows 10:14) correspond to g₁, which is a function
    % of the target position p_T.
    
    % Compute the derivative of g₁ = [τ₁; φ_Rx,₁; φ_Tx,₁; θ_Rx,₁; θ_Tx,₁] with respect to p_T.
    % Compute the two “look‐vectors” in the local coordinate systems:
    pT  = params.pT;
    pTx = params.pTx;
    pRx = params.pRx;
    
    % Tx and Rx local vectors:
    u1 = params.RTx' * (pT - pTx); % Tx–side vector.
    v1 = params.RRx' * (pT - pRx); % Rx–side vector.
    norm_u1 = norm(u1);
    norm_v1 = norm(v1);
    
    % Compute elevation angles from the local vectors:
    aodEl1 = asin(u1(3) / norm_u1);
    aoaEl1 = asin(v1(3) / norm_v1);
    
    % Standard basis vectors in ℝ³:
    e1 = [1; 0; 0];
    e2 = [0; 1; 0];
    e3 = [0; 0; 1];
    
    % Derivative of delay (τ₁) with respect to p_T (Eq. (70)):
    dtau_dpT = (1 / params.c) * ( (pT - pTx)' / norm(pT - pTx) + (pT - pRx)' / norm(pT - pRx) );
    
    % Derivative of the Rx–side azimuth (φ_Rx,₁) with respect to p_T (Eq. (71)):
    dphiRx_dpT = (v1(1)) / (v1(1)^2 + v1(2)^2) * e2.' * (eye(3) - (v1 * e1.') / v1(1)) * params.RRx.';
    
    % Derivative of the Tx–side azimuth (φ_Tx,₁) with respect to p_T (Eq. (72)):
    dphiTx_dpT = (u1(1)) / (u1(1)^2 + u1(2)^2) * e2.' * (eye(3) - (u1 * e1.') / u1(1)) * params.RTx.';
    
    % Derivative of the Rx–side elevation (θ_Rx,₁) with respect to p_T:
    daoaEl_dpT = (1 / (norm_v1 * cos(aoaEl1))) * e3.' * (eye(3) - (v1 * v1.') / (norm_v1^2)) * params.RRx.';
    
    % Derivative of the Tx–side elevation (θ_Tx,₁) with respect to p_T:
    daodEl_dpT = (1 / (norm_u1 * cos(aodEl1))) * e3.' * (eye(3) - (u1 * u1.') / (norm_u1^2)) * params.RTx.';
    
    % Assemble the 5×3 derivative block for g₁:
    % g₁ = [τ₁; φ_Rx,₁; φ_Tx,₁; θ_Rx,₁; θ_Tx,₁]
    % so that d(g₁)/d(p_T) is 5×3.
    dgdpt = [ dtau_dpT; dphiRx_dpT; dphiTx_dpT; daoaEl_dpT; daodEl_dpT ];  % 5×3
    
    % Place the derivative block into the lower block of the overall Jacobian.
    % Rows 10 to 14 of J correspond to g₁ and columns 10 to 12 correspond to p_T.
    J(10:14, 10:12) = dgdpt;

    % Add the one at (10, 3)
    J(10, 3) = 1;
end
