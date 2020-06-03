%% Adsorption rate constant
% T : Temperature (K)
% m : Mass of reactant (kg)

function k = k_ads(T, m)
    A = pi * (0.14 * 1e-9)^2;           % Site area (m2)
    kB = 1.38064852e-23;                % Boltzmann constant
    k = A / sqrt(2 * pi * m * kB * T);
end
