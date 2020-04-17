clc, clear, close all

%% Settings

% Plot defaults
set(0,'defaultAxesFontSize',12)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 10);

% Options
writeFiles = true;                 % Write files 
plotFigs = false;                  % Plot figures of transient behaviour
is_batch = true;                   % Batch or straight

% Free rate parameters
deactivation = 0.0;                % Amount of deactivated rxn with oxide
                                   %  - multiplies foward rxn rate
if (deactivation > 0)
    Ea_oxide = 136.5e3;            % Cu oxidation activation energy (J/mol)
else
    Ea_oxide = 204.4e3;            % Cu oxidation activation energy (J/mol)
end

% Parameters
T = 673;                           % Temperature (K)
P2 = 4;                            % System inlet pressure (bar)
pCO2 = 0;                          % CO2 inlet pressure (bar)
xO2s = 0:0.0002:0.005;             % O2 fraction added to Ar
xCO = 0.07;                        % CO fraction
atm = 101325;                      % Convert atm to Pa
tf = 3600;                         % Simulation time (s)

if ~is_batch
    % Straight
    ntanks = 10;                   % Number of tanks (CSTRs)
    if P2 == 2
        p = 1.2;                   % Pressure (bar) at 673 K
        F = 3.2e-13;               % Flow rate (m3/s) at 673 K
    else
        p = 2.5;                   % Pressure (bar) at 673 K
        F = 3.5e-13;               % Flow rate (m3/s) at 673 K
    end
    L = 216e-6;                    % Reactor length (m)
    D = 10e-6;                     % Reactor diameter (m)
    systID = 's';
else
    % Batch
    ntanks = 1;                    % Number of tanks (CSTRs)
    if P2 == 2
        p = 1.0;                   % Pressure (bar)
        F = 3.8e-13;               % Flow rate (m3/s)
    else
        p = 2.2;                   % Pressure (bar)
        F = 4.1e-13;               % Flow rate (m3/s)
    end
    L = 180e-6;                    % Reactor length (m)
    D = 120e-6;                    % Reactor diameter (m)
    systID = 'b';
end

% Reactor dimensions
h = 100e-9;                        % Height (m)
A = h * D;                         % Area (m2) pi * D^2
V = L * A;                         % Volume (m3)

% Catalyst features
rhoCat = 2e20;                     % Catalyst sites/m2
nCat = 1000;                       % Number of catalysts in arrays
cwidth = 120e-9;                   % Catalyst nanoparticle width (m)
ccollar = 40e-9 * pi * cwidth;     % Catalyst side (m2)
ccircle = pi * cwidth^2 * 0.25;    % Catalyst top (m2)
Acat = nCat * (ccircle + ccollar); % Total catalyst area (m2)

% Tanks in series model
Ltank = L / ntanks;                % Length of one tank (m)
tanks = 0:Ltank:L;                 % Incremental tank lengths (m)
Vtank = Ltank * A;                 % Volume of one tank (m3)
tau = Vtank / F;                   % Residence time (s)          
Acat = Acat / ntanks;              % Catalyst area in one tank (m2)

%% Run simulation(s)

disp('~~~~~~~~~~~~~~~~~~~~~')
disp('Simulation parameters')
disp('~~~~~~~~~~~~~~~~~~~~~')
disp(['Pressure = ' num2str(p) ' bar'])
disp(['Temperature = ' num2str(T) ' K'])
disp(['CO fraction = ' num2str(xCO)])
disp(['Residence time = ' num2str(tau) ' s'])
disp(['Final time = ' num2str(tf) ' s'])
disp(['Ntanks = ' num2str(ntanks)])
disp(['Ltank = ' num2str(Ltank) ' m'])
disp(['Acat = ' num2str(Acat) ' m2'])
disp(['Vtank = ' num2str(Vtank) ' m3'])
disp('~~~~~~~~~~~~~~~~~~~~~')
disp('Starting...')
disp('~~~~~~~~~~~~~~~~~~~~~')

% Initialise variables
yn = 7;                            % Number of variables
rn = 9;                            % Number of rate terms
pn = 3;                            % Number of gas phase species
y0 = zeros(yn, 1);                 % Solution vector (one tank)
y0ext = zeros(yn * ntanks, 1);
rates = zeros(length(xO2s), rn * ntanks);
concs = zeros(length(xO2s), pn * ntanks);
cover = zeros(length(xO2s), (yn - pn) * ntanks);
for j = 1:ntanks
    current = yn * (j - 1);
    y0ext(current + 1: current + yn) = y0;
end

% Loop over O2 fractions
for i = 1:length(xO2s)
    disp(['O2 frac: ' num2str(xO2s(i))])
    disp('---------------------')
    
    % Partial pressures in Pa
    ptot = p * 1e5;
    pCO = ptot * xCO;
    pO2 = ptot * xO2s(i);
    
    % System parameters required to solve
    params = [T; pCO; pO2; pCO2; rhoCat; F; Vtank; Acat; ntanks;...
        Ea_oxide; deactivation];
    
    fun = @(tvals, yvals)(get_CO_oxidation_odes(tvals, yvals, params));
    jac = @(tvals, yvals)(get_CO_oxidation_jac(tvals, yvals, params));
    
    % Set ODE options for ode15s
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'BDF', 'on',...
        'Stats', 'on', 'NonNegative', ones(yn * ntanks, 1));
    
    % Solve the ODEs
    tic
    [t, y] = ode15s(fun, [0, tf], y0ext, options);
    toc
    
    % Save final-time variables
    for j = 1 : ntanks
        current = yn * (j - 1);
        coverf = y(end, (current + 1):(current + yn - pn));
        concsf = y(end, (current + yn - pn + 1):(current + yn));
        
        % Coverages
        c_offset = (yn - pn) * (j - 1);
        cover(i, (c_offset + 1):(c_offset + yn - pn)) = coverf;
        
        % Pressures
        p_offset = pn * (j - 1);
        concs(i, (p_offset + 1):(p_offset + pn)) = concsf;
        
        % Rates
        r_offset = rn * (j - 1);
        rate_params = [T, Acat, Ea_oxide, deactivation];
        rates(i, (r_offset + 1):(r_offset + rn)) =...
            get_CO_oxidation_rates(coverf, rate_params); 
    end
    
    % Plot time dependent rates and coverages
    if plotFigs
        figure('units', 'normalized', 'outerposition', [0 0 1 1])
        set(gcf, 'color', 'white')
        for j = 1 : ntanks
            current = yn * (j - 1);
            allsites = 1 - sum(y(:, (current + 1):(current + yn - pn)), 2);
            subplot(2, ntanks, j)
            plot(t, y(:, current + 1), '-',...
                t, y(:, current + 2), '--',...
                t, y(:, current + 3), '-.',...
                t, y(:, current + 4), ':',...
                t, 1 - allsites, '-.')
            hold on
            ylabel('Coverage')
            set(gca, 'ylim', [0, 1])
            legend('CO', 'O_2', 'O', 'oxide', 'free', 'location', 'west')
            if ntanks > 1
                title(['xO_2 = ' num2str(xO2s(i)) ': tank ' num2str(j)])
            else
                title(['xO_2 = ' num2str(xO2s(i))])
            end

            subplot(2, ntanks, ntanks + j)
            plot(t, y(:, current + (yn - 2)), '-',...
                t, y(:, current + (yn - 1)), '--',...
                t, y(:, current + yn), '-.')
            xlabel('Time (s)')
            ylabel('Partial pressure (atm)')
            legend('CO', 'O_2', 'CO_2', 'location', 'west')
        end
    end

    disp('---------------------')
end

%% Results: Coverage vs CO2 concentration
figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
set(gcf, 'color', 'white')
hold on

% Plot CO2 partial pressure
yyaxis left
pO2s = concs(:, 2) * 1.01325;
pCO2s = concs(:, end) * 1.01325;
plot(xO2s * 100, pCO2s / p * 100, '-', 'color', [0 0 0.7])
if ntanks == 1
    plot(xO2s * 100, pO2s / p * 100, '-', 'color', [0.7 0 0])
end
ylabel('Concentration (%)', 'color', 'k')

% Plot coverage
yyaxis right
cls = colormap(copper);
step = round(size(cls, 1) / ntanks);
active_sites = 1 - cover(:, 4);
plot(xO2s * 100, active_sites, '-', 'color', cls(1, :))
for i = 2:ntanks
    i1 = (yn - 3) + (i - 1) * (yn - 3);
    i2 = step * (i - 1);
    active_sites = 1 - cover(:, i1);
    plot(xO2s * 100, active_sites, '-', 'color', cls(i2, :))
end
xlabel('O_2 concentration (%)')
ylabel('Active sites fraction')

% Format axes
h = gca;
h.YAxis(1).Color = 'k';
h.YAxis(2).Color = 'k';
h.YAxis(2).Limits = [0, 1];
set(gca,'XLim', [0 0.5])

% Add legend
child_handles = allchild(h);
if ntanks == 1
    l = legend([child_handles(end), child_handles(end-1),...
        child_handles(end-2)], 'CO_2', 'O_2', 'Coverage');
else
    l = legend([child_handles(end), child_handles(end-1)], 'CO_2',...
        'Coverage');
end
l.Box = 'Off';
l.Location = 'northoutside';
l.Orientation = 'horizontal';

%% File output
if writeFiles
    fending = [num2str(P2) 'bar_' num2str(ntanks) 'tanks_' systID];
    saveas(gcf, ['figs/cover_concs_' fending '.png'])
    csvwrite(['data/cover_' fending '.csv'], cover);
    csvwrite(['data/concs_' fending '.csv'], concs);
    csvwrite(['data/rates_' fending '.csv'], rates);
end

disp('Finished.')