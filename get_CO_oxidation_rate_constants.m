%% Elementary step reaction rates
% y      : Solution vector for surface coverage
%          y = (yCO, yO2, yO, yox)
% params : System parameters
%          params = (T, A, Ea_oxide, deactivation)

function rate_constants = get_CO_oxidation_rate_constants(params)
    % Parameters
    T = params(1);
    A = params(2);
    Ea_ox = params(3);              % Oxide activation energy (J/mol)
    deactivation = params(4);       % Deactivated reaction fraction

    % Constants
    mCO = 28 * 1.66054e-27;         % Molecular mass (a.u.)
    mO2 = 32 * 1.66054e-27;

    % Ads/des parameters from Filot (2018)
    sigmaCO = 1;                    % Symmetry number 
    sigmaO2 = 2;
    thetaCO = 2.8;                  % Rotational temperature (K)
    thetaO2 = 2.08;
    EadsCO = 80e3;                  % Desorption energy (J/mol)
    EadsO2 = 40e3;
    
    % Prefactors (transition state theory estimates)
    Afwd = 1e13;                    % (1/s)
    Arev = 1e13;

    % Activation energies from Falsig et al. (2008)
    Eafwd = 78.41e3;                % (J/mol)
    Eafwd_O = 12.94e3;
    Earev_O = 248.22e3;
        
    % Adsorption/desorption rate constants
    k_ads_CO = k_ads(T, A, mCO);
    k_ads_O2 = k_ads(T, A, mO2);    
    k_des_CO = k_des(T, A, mCO, sigmaCO, thetaCO, EadsCO);
    k_des_O2 = k_des(T, A, mO2, sigmaO2, thetaO2, EadsO2);

    % Reaction rate constants
    k_fwd = k_arr(T, Afwd, Eafwd);
    k_fwd_O = k_arr(T, Afwd, Eafwd_O);
    k_rev_O = k_arr(T, Arev, Earev_O);
    
    % Rate of oxide formation
    k_oxd = k_arr(T, Afwd, Ea_ox);
    
    % Rate of reaction with oxide
    k_fwd_ox = deactivation * k_fwd;
    
    % Combine to rate constant vector
    rate_constants = [k_ads_CO, k_ads_O2, k_des_CO, k_des_O2, k_fwd,...
        k_fwd_O, k_rev_O, k_oxd, k_fwd_ox]; 
end
