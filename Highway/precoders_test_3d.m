% precoders_test_3d.m

Nrx_x = 8;
Nrx_z = 8;
Ntx_x = 8;
Ntx_z = 8;

[W, F] = generate_precoders_3d(Ntx_x, Ntx_z, Nrx_x, Nrx_z, "QAMCOM");
[~, N] = size(W); 

figure(); 
imSize = 50;
rxPlot = imshow( ...
    zeros(imSize, imSize), ...
    "InitialMagnification", ...
    5000, ...
    'Colormap', ...
    parula);
clim([-20, 0])

figure();
txPlot = imshow( ...
    zeros(imSize, imSize), ...
    "InitialMagnification", ...
    5000, ...
    'Colormap', ...
    parula);
clim([-20, 0])

linkdata on

for n = 1:N

    w = reshape(W(:, n), [Nrx_z, Nrx_x]);
    rxPlot.CData = 10 * log10(fftshift(abs(fft2(w, imSize, imSize)) .^ 2));
    refreshdata
    drawnow limitrate

    f = reshape(F(:, n), [Ntx_z, Ntx_x]);
    txPlot.CData = 10 * log10(fftshift(abs(fft2(f, imSize, imSize)) .^ 2));
    refreshdata
    drawnow limitrate
end