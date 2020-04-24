clc, clear, close all

%% Settings

% Plot defaults
set(0, 'defaultAxesFontSize',12)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 10);

%% Optimise free parameters

% Initial guess
Ea_oxide = 100e3;                  % Cu oxidation activation energy (J/mol)
deactivation = 1.0;                % Amount of deactivated rxn with oxide
                                   %  - multiplies foward rxn rate

guess = [deactivation, Ea_oxide * 1e-5];
disp('Starting optimisation with initial values:')
disp(['* Ea_oxide: ' num2str(Ea_oxide * 1e-3, '%1.0f') ' kJ/mol'])
disp(['* deactivation: ' num2str(deactivation, '%1.2f')])
disp('------------------------------------')

% Optimisation
options = optimoptions('fsolve', 'Display', 'iter', 'UseParallel', true);
[x, fval, flag, output] = fsolve(@parameter_optimisation, guess, options);
disp('Optimisation output:')
disp('* Values:')
disp(fval)
disp('* Flag:')
disp(flag)
disp('------------------------------------')

% Final values
Ea_oxide = abs(x(2)) * 1e5;
deactivation = abs(x(1));
disp('Final values:')
disp(['* Ea_oxide: ' num2str(Ea_oxide * 1e-3, '%1.0f') ' kJ/mol'])
disp(['* deactivation: ' num2str(deactivation, '%1.2f')])
disp('------------------------------------')

% options = optimset('Display', 'iter', 'UseParallel', true,...
%     'PlotFcns', @optimplotfval);
% fminsearch(@parameter_optimisation, guess, options)