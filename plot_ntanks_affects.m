clc, clear, close all

%% Settings

% Plot defaults
set(0,'defaultAxesFontSize',12)
set(0, 'DefaultLineLineWidth', 1);
set(0, 'DefaultLineMarkerSize', 10);

xO2s = 0:0.0002:0.005;
xO2s = xO2s * 100;

figure('PaperUnits', 'inches', 'PaperPosition', [0 0 5 3.3])
set(gcf, 'color', 'white')
hold on

figure('PaperUnits', 'inches', 'PaperPosition', [0 0 10 5])
set(gcf, 'color', 'white')
hold on
cls = colormap(copper);
step = round(size(cls, 1) / 10);

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