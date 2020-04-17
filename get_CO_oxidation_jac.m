%% Jacobian for tanks-in-series model with surface kinetics
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

function j = get_CO_oxidation_jac(~, y, params)

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

    % Initialise jacobian matrix
    ny = length(y) / n;
    j = zeros(n * ny, n * ny);

    % Loop over tanks
    for i = 1:n
        % Variables for current tank
        current = ny * (i - 1);        % Offset to current tank indices
        previous = ny * (i - 2);       % Offset to previous tank indices
        yCO = y(current + 1);          % CO coverage
        yO2 = y(current + 2);          % O2 coverage
        yO = y(current + 3);           % O coverage
        yox = y(current + 4);          % ox coverage
        pCO = y(current + (ny - 2));   % pCO (atm)
        pO2 = y(current + (ny - 1));   % pO2 (atm)
        pCO2 = y(current + ny);        % pCO2 (atm)
        
        % Free sites
        yf = 1 - sum(y(current + 1 : current + (ny - 3)));
        
        % Get process rate constants
        rate_params = [T, A, Ea_oxide, deactivation];
        rate_constants = get_CO_oxidation_rate_constants(rate_params);
        k_ads_CO = rate_constants(1);
        k_ads_O2 = rate_constants(2);
        k_des_CO = rate_constants(3);
        k_des_O2 = rate_constants(4);
        k_fwd = rate_constants(5);
        k_fwd_O = rate_constants(6);
        k_rev_O = rate_constants(7);
        k_oxd = rate_constants(8);
        k_fwd_ox = rate_constants(9);
        
        % Row/column indices for variables in each tank:
        % 1-7, 8-14, 15-21, 22-28, 29-35, 36-42, 43-49, 50-56, 57-63, 64-70
        
        % ddyCO/d...
        j(1 + current, 1 + current) = -k_ads_CO * pCO * atm - k_des_CO -...
            k_fwd * yO - k_fwd_ox * yox;
        j(1 + current, 2 + current) = -k_ads_CO * pCO * atm;
        j(1 + current, 3 + current) = -k_ads_CO * pCO * atm - k_fwd * yCO;
        j(1 + current, 4 + current) = -k_ads_CO * pCO * atm -...
            k_fwd_ox * yCO;
        j(1 + current, 5 + current) = k_ads_CO * yf;
        j(1 + current, 6 + current) = 0;
        j(1 + current, 7 + current) = 0;
        
        % ddyO2/d...
        j(2 + current, 1 + current) = -k_ads_O2 * pO2 * atm +...
            k_fwd_O * yO2;
        j(2 + current, 2 + current) = -k_ads_O2 * pO2 * atm - k_des_O2 -...
            k_fwd_O * yf + k_fwd_O * yO2;
        j(2 + current, 3 + current) = -k_ads_O2 * pO2 * atm +...
            k_fwd_O * yO2 + 2 * k_rev_O * yO;
        j(2 + current, 4 + current) = -k_ads_O2 * pO2 * atm +...
            k_fwd_O * yO2;
        j(2 + current, 5 + current) = 0;
        j(2 + current, 6 + current) = k_ads_O2 * yf;
        j(2 + current, 7 + current) = 0;
        
        % ddyO/d...
        j(3 + current, 1 + current) = -2 * k_fwd_O * yO2 - k_fwd * yO;
        j(3 + current, 2 + current) = -2 * k_fwd_O * yO2 +...
            2 * k_fwd_O * yf;
        j(3 + current, 3 + current) = -2 * k_fwd_O * yO2 -...
             4 * k_rev_O * yO - k_fwd * yCO - k_oxd;
        j(3 + current, 4 + current) = -2 * k_fwd_O * yO2;
        j(3 + current, 5 + current) = 0;
        j(3 + current, 6 + current) = 0;
        j(3 + current, 7 + current) = 0;
        
        % ddyox/d...
        j(4 + current, 1 + current) = - k_fwd_ox * yox;  
        j(4 + current, 2 + current) = 0;
        j(4 + current, 3 + current) = k_oxd;
        j(4 + current, 4 + current) = - k_fwd_ox * yCO;
        j(4 + current, 5 + current) = 0;
        j(4 + current, 6 + current) = 0;
        j(4 + current, 7 + current) = 0;
        
        % ddpCO/d...
        j(5 + current, 1 + current) = sigma * k_des_CO + sigma *...
            k_ads_CO * pCO * atm;
        j(5 + current, 2 + current) = sigma * k_ads_CO * pCO * atm;
        j(5 + current, 3 + current) = sigma * k_ads_CO * pCO * atm;
        j(5 + current, 4 + current) = sigma * k_ads_CO * pCO * atm;
        j(5 + current, 5 + current) = -itau * atm - sigma * k_ads_CO * yf;
        j(5 + current, 6 + current) = 0;
        j(5 + current, 7 + current) = 0;
        
        % ddpO2/d...
        j(6 + current, 1 + current) = sigma * k_ads_O2 * pO2 * atm;
        j(6 + current, 2 + current) = sigma * k_des_O2 + sigma *...
            k_ads_O2 * pO2 * atm;
        j(6 + current, 3 + current) = sigma * k_ads_O2 * pO2 * atm;
        j(6 + current, 4 + current) = sigma * k_ads_O2 * pO2 * atm;
        j(6 + current, 5 + current) = 0;
        j(6 + current, 6 + current) = -itau * atm - sigma * k_ads_O2 * yf;
        j(6 + current, 7 + current) = 0;
        
        % ddpCO2/d...
        j(7 + current, 1 + current) = sigma * k_fwd * yO + sigma *...
            k_fwd_ox * yox;
        j(7 + current, 2 + current) = 0;
        j(7 + current, 3 + current) = sigma * k_fwd * yCO;
        j(7 + current, 4 + current) = sigma * k_fwd_ox * yCO;
        j(7 + current, 5 + current) = 0;
        j(7 + current, 6 + current) = 0;
        j(7 + current, 7 + current) = -itau * atm;
        
        % Add terms for inlet pressure (outlet from previous tank)
        if i > 1
            j(5 + current : 7 + current, 5 + previous : 7 + previous) =...
                diag(ones(3, 1) * itau);
        end
       
       % Scale to match odes for pressure in atm
       j(5 + current : 7 + current, 1 + current : 7 + current) = ...
           j(5 + current : 7 + current, 1 + current : 7 + current) / atm;   
       j(1 + current : 7 + current, 5 + current : 7 + current) = ...
           j(1 + current : 7 + current, 5 + current : 7 + current) * atm; 
    end
end