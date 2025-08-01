% pattern_test.m

phis = linspace(-pi / 2, pi / 2, 200);
thetas = linspace(-pi / 2, pi / 2, 200);

pattern = zeros(200, 200);
for i = 1:200
    for j = 1:200
        pattern(i, j) = antenna_element_pattern(phis(i), thetas(j), 65 / 180 * pi, 8);
    end
end

figure()
imshow(10 * log10(pattern .^2), "Colormap", parula, "InitialMagnification", 1000);
clim([-20, 8]);
colorbar