clc, clear, close all

%% Settings

% Plot defaults
set(0,'defaultAxesFontSize',12)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 10);

% Options 
writeFiles = false;                % Write files 
plotFigs = true;                   % Plot figures of transient behaviour

% Free rate parameters
deactivation = 0.2;                % Amount of deactivated rxn with oxide
                                   %  - multiplies foward rxn rate
if (deactivation > 0)
    Ea_oxide = 136.5e3;            % Cu oxidation activation energy (J/mol)
else
    Ea_oxide = 204.4e3;            % Cu oxidation activation energy (J/mol)
end

% Parameters
T = 673;                           % Temperature (K)
p = 4;                             % Pressure (bar)
pCO2 = 0;                          % CO2 inlet pressure (bar)
xO2s = 0:0.001:0.005;              % O2 fraction added to Ar
xCO = 0.07;                        % CO fraction
atm = 101325;                      % Convert atm to Pa
tf = 3600;                         % Simulation time (s)

% Catalyst features
nCat = 1;                          % Number of catalysts in arrays
cwidth = 120e-9;                   % Catalyst nanoparticle width (m)
ccollar = 40e-9 * pi * cwidth;     % Catalyst side (m2)
ccircle = pi * cwidth^2 * 0.25;    % Catalyst top (m2)
Acat = nCat * (ccircle + ccollar); % Total catalyst area (m2)

%% Run simulation(s)

% Display parameters
disp('~~~~~~~~~~~~~~~~~~~~~')
disp('Simulation parameters')
disp('~~~~~~~~~~~~~~~~~~~~~')
disp(['Pressure = ' num2str(p) ' bar'])
disp(['Temperature = ' num2str(T) ' K'])
disp(['CO fraction = ' num2str(xCO)])
disp(['Final time = ' num2str(tf) ' s'])
disp(['Acat = ' num2str(Acat) ' m2'])
disp('~~~~~~~~~~~~~~~~~~~~~')
disp('Starting...')
disp('~~~~~~~~~~~~~~~~~~~~~')

% Initialise variables
yn = 4;                            % Number of variables
rn = 9;                            % Number of rate terms
y0 = zeros(yn, 1);                 % Solution vector
rates = zeros(length(xO2s), rn);
cover = zeros(length(xO2s), yn);

% Loop over O2 fractions
for i = 1:length(xO2s)
    disp(['O2 frac: ' num2str(xO2s(i))])
    disp('---------------------')
    
    % Partial pressures in Pa
    ptot = p * 1e5;
    pCO = ptot * xCO;
    pO2 = ptot * (1 - xCO) * xO2s(i);
    
    % System parameters required to solve
    params = [T; pCO; pO2; pCO2; Acat; Ea_oxide; deactivation];
    
    fun = @(tvar, yvar)(get_CO_oxidation_surface_odes(tvar, yvar, params));
    
    % Set ODE options for ode15s
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'BDF', 'on',...
        'Stats', 'on', 'NonNegative', ones(yn, 1)); 
    
    % Solve the ODEs
    tic
    [t, y] = ode15s(fun, [0, tf], y0, options);
    toc
    
    % Save final-time variables      
    cover(i, :) = y(end, :);
    rate_params = [T, Acat, Ea_oxide, deactivation];
    rates(i, :) = get_CO_oxidation_rates(y(end, :), rate_params); 
        
    % Plot time-dependent rates and coverages
    if plotFigs
        figure('units', 'normalized', 'outerposition', [0 0 3.3 3.3])
        set(gcf, 'color', 'white')
        plot(t, y(:, 1), '-', t, y(:, 2), '--', t, y(:, 3), '-.',...
            t, y(:, 4), ':', t, 1 - sum(y, 2), '-.')
        ylabel('Coverage')
        set(gca, 'ylim', [0, 1])
        legend('CO', 'O_2', 'O', 'oxide', 'free', 'location', 'west')
        title(['xO_2 = ' num2str(xO2s(i))])
    end

    disp('---------------------')
end

%% Results: Coverage vs CO2 turnover rate
figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
set(gcf, 'color', 'white')
hold on

% Plot CO2 production rate
yyaxis left
rCO2s = rates(:, 5) + rates(:, 9);
plot(xO2s * 100, rCO2s, '-x', 'markersize', 6, 'color', [0 0 0.7])
ylabel('TOF (1/(sites.s))', 'color', 'k')

% Plot coverage
yyaxis right
active_sites = 1 - cover(:, end);
plot(xO2s * 100, active_sites, '-+', 'color', [0.7 0 0])
xlabel('O_2 concentration (%)')
ylabel('Active sites fraction')

% Format axes
h = gca;
h.YAxis(1).Color = 'k';
h.YAxis(2).Color = 'k';
h.YAxis(2).Limits = [0, 1];

% Add legend
child_handles = allchild(h);
l = legend([child_handles(end), child_handles(1)], 'TOF', 'Sites');
l.Box = 'Off';
l.Location = 'northoutside';
l.Orientation = 'horizontal';

%% File output
if writeFiles
    fending = [num2str(P2) 'bar_' num2str(T) 'K'];
    saveas(gcf, ['figs/cover_rate_' fending '.png'])
    csvwrite(['data/cover_' fending '.csv'], cover);
    csvwrite(['data/concs_' fending '.csv'], concs);
    csvwrite(['data/rates_' fending '.csv'], rates);
end

disp('Finished.')