%% Parameter optimisation for CO oxidation
% opt_vec : Optimisable parameters
%           opt_vec = (deactivation, Ea_oxide)

function f = parameter_optimisation(opt_vec)

f = zeros(2, 1);

% Free rate parameters
deactivation = abs(opt_vec(1));    % Amount of deactivated rxn with oxide
                                   %  - multiplies foward rxn rate
Ea_oxide = abs(opt_vec(2)) * 1e5;  % Cu oxidation activation energy (J/mol)

% Parameters
T = 673;                           % Temperature (K)
P2 = 4;                            % System inlet pressure (bar)
pCO2 = 0;                          % CO2 inlet pressure (bar)
xO2s = 0:0.0001:0.005;             % O2 fraction added to Ar
xCO = 0.07;                        % CO fraction
tf = 3600;                         % Simulation time (s)

% Batch
ntanks = 1;                        % Number of tanks (CSTRs)
if P2 == 2
    p = 1.0;                       % Pressure (bar)
    F = 7.3e-13;                   % Flow rate (m3/s)
else
    p = 2.2;                       % Pressure (bar)
    F = 7.5e-13;                   % Flow rate (m3/s)
end
L = 180e-6;                        % Reactor length (m)
D = 120e-6;                        % Reactor diameter (m)

% Reactor dimensions
h = 100e-9;                        % Height (m)
A = h * D;                         % Area (m2) pi * D^2

% Catalyst features
rhoCat = 6e19;                     % Catalyst sites/m2
nCat = 1000;                       % Number of catalysts in arrays
cwidth = 120e-9;                   % Catalyst nanoparticle width (m)
ccollar = 40e-9 * pi * cwidth;     % Catalyst side (m2)
ccircle = pi * cwidth^2 * 0.25;    % Catalyst top (m2)
Acat = nCat * (ccircle + ccollar); % Total catalyst area (m2)

% Tanks in series model
Ltank = L / ntanks;                % Length of one tank (m)
Vtank = Ltank * A;                 % Volume of one tank (m3)
Acat = Acat / ntanks;              % Catalyst area in one tank (m2)

%% Run simulation(s)

% Initialise variables
yn = 7;                            % Number of variables
pn = 3;                            % Number of gas phase species
y0 = zeros(yn, 1);                 % Solution vector (one tank)
y0ext = zeros(yn * ntanks, 1);
cover = zeros(length(xO2s), (yn - pn) * ntanks);
concs = zeros(length(xO2s), pn * ntanks);
for j = 1:ntanks
    current = yn * (j - 1);
    y0ext(current + 1: current + yn) = y0;
end

% Loop over O2 fractions
for i = 1:length(xO2s)
    
    % Partial pressures in Pa
    ptot = p * 1e5;
    pCO = ptot * xCO;
    pO2 = ptot * xO2s(i);
    
    % System parameters required to solve
    params = [T; pCO; pO2; pCO2; rhoCat; F; Vtank; Acat; ntanks;...
        Ea_oxide; deactivation];
    
    fun = @(tvals, yvals)(get_CO_oxidation_odes(tvals, yvals, params));
    
    % Set ODE options for ode15s
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'BDF', 'on',...
        'NonNegative', ones(yn * ntanks, 1));
    
    % Solve the ODEs
    [~, y] = ode15s(fun, [0, tf], y0ext, options);
    
    % Save final-time variables
    for j = 1 : ntanks
        current = yn * (j - 1);
        coverf = y(end, (current + 1):(current + yn - pn));
        concsf = y(end, (current + yn - pn + 1):(current + yn));
        
        c_offset = (yn - pn) * (j - 1);
        cover(i, (c_offset + 1):(c_offset + yn - pn)) = coverf;
        
        p_offset = pn * (j - 1);
        concs(i, (p_offset + 1):(p_offset + pn)) = concsf;
    end
end

target1 = 0.003;
imax = find(concs(:, end) == max(concs(:, end)));
xO2max = xO2s(imax);
f(1) = (xO2max - target1) * 100;

target2 = 0.2;
xCO2max = concs(imax, end) * 1.01325 / p;
xCO2fin = concs(end, end) * 1.01325 / p;
Xmax = 0.5 * xCO2max / xO2max;
Xfin = 0.5 * xCO2fin / xO2s(end);
f(2) = (Xfin / Xmax - target2);
 
% target3 = 0;
% active_sites = sum(cover(end, end - 3:end - 1));
% f(3) = abs(active_sites - target3);


% target1 = 0.003;
% imax = find(concs(:, end) == max(concs(:, end)));
% xO2max = xO2s(imax);
% fa = (xO2max - target1)^2 * 100;
% 
% target2 = 0.16;
% xCO2max = concs(imax, end) * 1.01325 / p;
% xCO2fin = concs(end, end) * 1.01325 / p;
% Xmax = 0.5 * xCO2max / xO2max;
% Xfin = 0.5 * xCO2fin / xO2s(end);
% fb = (Xfin - target2)^2;
%  
% target3 = 0.8;
% active_sites = sum(cover(end, end - 3:end - 1));
% fc = (Xmax - target3)^2;
% 
% fd = active_sites;
% 
% f = fa + fb + fc;


end