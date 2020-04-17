%% Elementary step reaction rates
% y      : Solution vector for surface coverage
%          y = (yCO, yO2, yO, yox)
% params : System parameters
%          params = (T, A, Ea_oxide, deactivation)

% Reactions
% =========
% Overall: 2CO + O2 -> 2CO2    (f/r)
% 1) CO + *   <-> CO*          (1/2)
% 2) O2 + *   <-> O2*          (3/4)
% 3) O2* + *  <-> 2O*          (5/6)
% 4) CO* + O*  -> CO2 + 2*     (7)
% 5) O*        -> Oox          (8)
% 6) CO* + Oox -> CO2 + *      (9)

function rates = get_CO_oxidation_rates(y, params)
    % Get rate constants
    rate_constants = get_CO_oxidation_rate_constants(params);
    k_ads_CO = rate_constants(1);
    k_ads_O2 = rate_constants(2);
    k_des_CO = rate_constants(3);
    k_des_O2 = rate_constants(4);
    k_fwd = rate_constants(5);
    k_fwd_O = rate_constants(6);
    k_rev_O = rate_constants(7);
    k_oxd = rate_constants(8);
    k_fwd_ox = rate_constants(9);

    % Free sites
    yfree = 1 - sum(y);
    
    % Adsorption/desorption rates
    adsCO = k_ads_CO * yfree;
    desCO = k_des_CO * y(1);
    adsO2 = k_ads_O2 * yfree;
    desO2 = k_des_O2 * y(2);

    % Reaction rate
    fwdRN = k_fwd * y(1) * y(3);
    
    % Dissociation/recombination rates
    fwdRN_O = k_fwd_O * y(2) * yfree;
    revRN_O = k_rev_O * y(3) * y(3);

    % Oxide formation rate
    OtoOox = k_oxd * y(3);
    
    % Reaction with oxide
    fwdRNox = k_fwd_ox * y(1) * y(4);

    % Combine to rates vector
    rates = [adsCO, adsO2, desCO, desO2, fwdRN, fwdRN_O, revRN_O,...
        OtoOox, fwdRNox];
end