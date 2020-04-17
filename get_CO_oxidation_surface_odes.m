%% Surface dynamics ODEs
% t      : Time (s), not used but required for ode solver
% y      : Solution vector for surface coverage
%          y = (yCO, yO2, yO, yox)
% params : System parameters
%          params = (T, pCO, pO2, pCO2, A, Ea_oxide, deactivation)

% Reactions
% =========
% Overall: 2CO + O2 -> 2CO2
% 1) CO + *   <-> CO*
% 2) O2 + *   <-> O2*
% 3) O2* + *  <-> 2O*
% 4) CO* + O*  -> CO2 + 2*
% 5) O*        -> Oox
% 6) CO* + Oox -> CO2 + *

function dy_dt = get_CO_oxidation_surface_odes(~, y, params)
    % Parameters
    T = params(1);
    pCO = params(2);
    pO2 = params(3);
    A = params(5);
    Ea_oxide = params(6);
    deactivation = params(7);
    
    % Get process rates
    rate_params = [T, A, Ea_oxide, deactivation];
    rates = get_CO_oxidation_rates(y, rate_params);
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

    % Combine to return ODE vector
    dy_dt = [dyCO_dt; dyO2_dt; dyO_dt; dyox_dt];
end