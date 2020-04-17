%% Adsorption rate constant
% T  : Temperature (K)
% P  : Pressure (Pa)
% SA : Surface area (m^2)
% m  : Mass of reactant (kg)

function k = k_ads(T, SA, m) % P, 
    kB = 1.38064852e-23;            % Boltzmann constant
    k = (SA) / sqrt(2 * pi * m * kB * T); %P * 
end