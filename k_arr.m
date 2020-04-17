%% Surface reaction rate constant
% T  : Temperature (K)
% A  : Pre-exponential factor (1/s)
% Ea : Activation energy (J/mol)

function k = k_arr(T, A, Ea)
    R = 8.3144598;                  % Gas constant
    k = A * exp(-Ea / (R * T));
end