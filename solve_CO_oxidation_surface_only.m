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
deactivation = 1.33;               % Amount of deactivated rxn with oxide
                                   %  - multiplies foward rxn rate
Ea_oxide = 84e3;                   % Cu oxidation activation energy (J/mol)

% Parameters
T = 673;                           % Temperature (K)
p = 2;                             % Pressure (bar)
pCO2 = 0;                          % CO2 inlet pressure (bar)
xO2s = 0:0.001:0.005;              % O2 fraction added to Ar
xCO = 0.07;                        % CO fraction
atm = 101325;                      % Convert atm to Pa
tf = 0.05;                         % Simulation time (s)

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
% y0(4) = 0.8;                     % Change to start with oxidised catalyst
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
        figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
        set(gcf, 'color', 'white')
        hold on
        plot(t * 1e3, y(:, 1), '-', 'color', 'r')
        plot(t * 1e3, y(:, 2), '--', 'color', [0.7 0.7 0])
        plot(t * 1e3, y(:, 3), '--', 'color', [0 0 0.7])
        plot(t * 1e3, y(:, 4), '-.', 'color', [0 0.7 0])
        plot(t * 1e3, 1 - sum(y, 2), 'k:')
        ylabel('Fraction')
        xlabel('Time (ms)')
        set(gca, 'ylim', [0, 1])
        l = legend('CO', 'O_2', 'O', 'Oxide', 'Free', 'location', 'east');  
        l.Box = 'Off';
        box on
        set(gca, 'LineWidth', 2)
        title(['xO_2 = ' num2str(xO2s(i))])
        saveas(gcf, 'figs/All_cover.png')
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
plot(xO2s * 100, rCO2s, '-', 'color', [0 0 0.7])
ylabel('TOF (1/(sites.s))', 'color', 'k')

% Plot coverage
yyaxis right
active_sites = 1 - cover(:, end);
plot(xO2s * 100, active_sites, '-', 'color', [0.7 0 0])
xlabel('O_2 concentration (%)')
ylabel('Active fraction')

% Format axes
h = gca;
h.YAxis(1).Color = 'k';
h.YAxis(2).Color = 'k';
h.YAxis(2).Limits = [0, 1];

% Add legend
child_handles = allchild(h);
l = legend([child_handles(end), child_handles(1)], 'TOF', 'Sites');
l.Box = 'Off';
l.Location = 'north';
l.Orientation = 'horizontal';

box on
h.LineWidth = 2;

%% File output
if writeFiles
    fending = [num2str(p) 'bar_' num2str(T) 'K'];
    saveas(gcf, ['figs/cover_rate_' fending '.png'])
    csvwrite(['data/cover_' fending '.csv'], cover);
    csvwrite(['data/concs_' fending '.csv'], concs);
    csvwrite(['data/rates_' fending '.csv'], rates);
end

disp('Finished.')