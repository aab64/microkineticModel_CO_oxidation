%% Desorption rate constant
% T      : Temperature (K)
% m      : Mass of reactant (kg)
% sigma  : Symmetry number
% theta  : Rotational temperature (K)
% Edes   : Desorption energy (J/mol)

function k = k_des(T, m, sigma, theta, Edes)
    A = pi * (0.14 * 1e-9)^2;       % Site area (m2)
    R = 8.3144598;                  % Gas constant
    kB = 1.38064852e-23;            % Boltzmann constant
    h = 6.62607004e-34;             % Planck constant
    k = kB * (T^3) / (h^3) * A * (2 * pi * m * kB) / ...
        (sigma * theta) * exp(-Edes / (R * T));
end