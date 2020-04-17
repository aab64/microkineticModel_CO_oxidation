%% Adsorption rate constant
% T : Temperature (K)
% A : Surface area (m^2)
% m : Mass of reactant (kg)

function k = k_ads(T, A, m)
    kB = 1.38064852e-23;               % Boltzmann constant
    k = A / sqrt(2 * pi * m * kB * T); 
end
