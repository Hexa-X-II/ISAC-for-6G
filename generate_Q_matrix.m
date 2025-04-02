function Q = generate_Q_matrix(Nx, Nz, lambda)
% Q = generate_Q_matrix(Nx, Nz, lambda)
%   Generate the Q matrix in alignment with precoders for the 3D case.

xs = -(Nx - 1) / 2:(Nx - 1) / 2;
zs = -(Nz - 1) / 2:(Nz - 1) / 2;

Q = zeros(3, Nx * Nz);
count = 1;
for x = xs
    for z = zs
        Q(:, count) = lambda / 2 * [x; 0; z];
        count = count + 1;
    end
end

end