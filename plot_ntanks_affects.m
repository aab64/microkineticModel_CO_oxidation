clc, clear, close all

%% Settings

% Plot defaults
set(0,'defaultAxesFontSize',12)
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultLineMarkerSize', 10);

xO2s = 0:0.0001:0.005;
xO2s = xO2s * 100;

figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
set(gcf, 'color', 'white')
hold on

figure('PaperUnits', 'inches', 'PaperPosition', [0 0 10 5])
set(gcf, 'color', 'white')
hold on
cls = colormap(copper);
step = round(size(cls, 1) / 10);

xO2maxs = zeros(10, 1);
Xmaxs = zeros(10, 1);

for ntanks = 1:10
    fending = ['4bar_' num2str(ntanks) 'tanks_s'];
    cover = csvread(['data/oxide_rxn/cover_' fending '.csv']);
    concs = csvread(['data/oxide_rxn/concs_' fending '.csv']);
    
    pCOs = concs(:, end - 2) * 101.325 / 2.5;
    pO2s = concs(:, end - 1) * 101.325 / 2.5;
    pCO2s = concs(:, end) * 101.325 / 2.5;   
    
    figure(1)
    plot(xO2s, pCO2s, 'color', cls(1 + step * (ntanks - 1), :)) 

    figure(2)
    subplot(2, 5, ntanks)
    hold on
    active_sites = 1 - cover(:, 4);
    plot(xO2s, active_sites, 'color', cls(1, :))
    for i = 2:ntanks
        i1 = 4 + (i - 1) * 4;
        i2 = step * (i - 1);
        active_sites = 1 - cover(:, i1);
        plot(xO2s, active_sites, 'color', cls(i2, :))
    end
    if (ntanks > 5)
        xlabel('O_2 concentration (%)')
    end
    if (ntanks == 1) || (ntanks == 6)
        ylabel('Active site fraction')
    end
    title(['n=' num2str(ntanks)])
    set(gca, 'xlim', [0 0.3])
    
    disp('=====================')
    imax = find(concs(:, end) == max(concs(:, end)));
    disp(['Onset at xO2 = ' num2str(xO2s(imax) * 100, '%1.3f') '%']);
    
    xO2maxs(ntanks) = xO2s(imax) * 100;
    
    xCO2max = concs(imax, end) * 1.01325 / 2.5;
    xCO2fin = concs(end, end) * 1.01325 / 2.5;
    Xmax = 0.5 * xCO2max / xO2s(imax);
    Xfin = 0.5 * xCO2fin / xO2s(end);
    disp(['Final conversion = ' num2str(100 * Xfin / Xmax, '%1.3f') '%']);
    
    Xmaxs(ntanks) = 100 * Xfin / Xmax;
    
    active_sites = 1 - cover(end, end);
    disp(['Final active site fraction = ' num2str(active_sites, '%1.3f')])
    disp('=====================')
end

figure(1)
xlabel('O_2 concentration (%)')
ylabel('CO_2 concentration (%)')
l = legend('1','2','3','4','5','6','7','8','9','10');
l.Box = 'off';
l.Location = 'eastoutside';
t = text(0.5, 0.68, 'Number of tanks', 'fontsize', 10);
saveas(gcf, ['figs/co2_conc_vs_ntanks_' fending '.png'])

figure(2)
saveas(gcf, ['figs/cover_vs_ntanks_' fending '.png'])

figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
set(gcf, 'color', 'white')
hold on
yyaxis left
plot(1:10, xO2maxs, '-o', 'color', [0.7 0 0], 'Linewidth', 2)
ylabel('O_2 inlet with max conversion (%)')
yyaxis right
plot(1:10, Xmaxs, '-^', 'color', [0 0 0.7], 'Linewidth', 2)
ylabel('Conversion at max O_2 inlet (%)')
xlabel('Number of tanks')
h = gca;
h.YAxis(1).Color = 'k';
h.YAxis(2).Color = 'k';
q1 = quiver(2.2, 23.5, -1, 0, 0, 'Linewidth', 2, 'color', [0.7 0 0], 'MaxHeadSize', 2);
q1.LineStyle = '-';
q2 = quiver(7, 24, 1, 0, 0, 'Linewidth', 2, 'color', [0 0 0.7], 'MaxHeadSize', 2);
q2.LineStyle = '-';
set(gca, 'YTick', 23:27)
saveas(gcf, ['figs/params_vs_ntanks_' fending '.png'])
