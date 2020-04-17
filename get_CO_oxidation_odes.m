%% CSTR rate equations and surface dynamics
% t      : Time (s), not used but required for ode solver
% y      : Solution vector for surface coverage and pressure
%          y_tanki = (yCO, yO2, yO, yox, pO2, pCO, pCO2)
%          y = zeros(n * 7, 1);
% params : System parameters
%          params = (T, pCO, pO2, pCO2, rhoCat, F, Vtank, A, ntanks,...
%                    Ea_oxide, deactivation)

% Reactions
% =========
% Overall: 2CO + O2 -> 2CO2
% 1) CO + *   <-> CO*
% 2) O2 + *   <-> O2*
% 3) O2* + *  <-> 2O*
% 4) CO* + O*  -> CO2 + 2*
% 5) O*        -> Oox
% 6) CO* + Oox -> CO2 + *

function dy_dt = get_CO_oxidation_odes(~, y, params)
    % Parameters
    T = params(1);
    pCOin = params(2);
    pO2in = params(3);
    pCO2in = params(4);
    rhoCat = params(5);
    F = params(6);
    V = params(7);
    A = params(8);
    n = params(9);
    Ea_oxide = params(10);
    deactivation = params(11);

    % Constants
    atm = 101325;                      % Convert atm to Pa
    kB = 1.38064852e-23;               % Boltzmann constant
    sigma = kB * T * rhoCat * A / V;   % Convert rate to Pa/s
    itau = F / V;                      % Inverse residence time (1/s)                     

    % Initialise ode vector
    ny = length(y) / n;
    dy_dt = zeros(n * ny, 1);

    % Loop over tanks in series
    for i = 1:n
        % Variables for current tank
        current = ny * (i - 1);        % Offset to current tank indices
        previous = ny * (i - 2);       % Offset to previous tank indices
        if i > 1
            pCOin = y(previous + (ny - 2)) * atm;
            pO2in = y(previous + (ny - 1)) * atm;
            pCO2in = y(previous + ny) * atm;
        end
        pCO = y(current + (ny - 2)) * atm;
        pO2 = y(current + (ny - 1)) * atm;
        pCO2 = y(current + ny) * atm;
        covers = y(current + 1 : current + (ny - 3));
        
        % Get process rates
        rate_params = [T, A, Ea_oxide, deactivation];
        rates = get_CO_oxidation_rates(covers, rate_params);
        adsCO = rates(1) * pCO;
        adsO2 = rates(2) * pO2;
        desCO = rates(3);
        desO2 = rates(4);
        fwdRN = rates(5);
        fwdRN_O = rates(6);
        revRN_O = rates(7);
        OtoOox = rates(8);
        fwdRNox = rates(9);

        % Coverage ODEs
        dyCO_dt = adsCO - desCO - fwdRN - fwdRNox;
        dyO2_dt = adsO2 - desO2 - fwdRN_O + revRN_O;
        dyO_dt = 2 * fwdRN_O - 2 * revRN_O - fwdRN - OtoOox; 
        dyox_dt = OtoOox - fwdRNox;
        
        % Pressure ODEs
        dpCO_dt = (pCOin - pCO) * itau + sigma * (desCO - adsCO); 
        dpO2_dt = (pO2in - pO2) * itau + sigma * (desO2 - adsO2);
        dpCO2_dt = (pCO2in - pCO2) * itau + sigma * (fwdRN + fwdRNox);

        % Rescale to match coverage magnitude for stability
        dpCO_dt = dpCO_dt / atm;
        dpO2_dt = dpO2_dt / atm;
        dpCO2_dt = dpCO2_dt / atm;

        % Combine to return ODE vector
        dy_dt(current + 1: current + ny) = [dyCO_dt; dyO2_dt;...
            dyO_dt; dyox_dt; dpCO_dt; dpO2_dt; dpCO2_dt];
    end
end