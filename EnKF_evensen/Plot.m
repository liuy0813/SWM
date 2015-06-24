
%% Plot results

% Plot size
FST = 40;           % Font size for Title
FS = 20;            % Font size for axes
LF = 20;            % Legend size
MS = 10;            % Marker size
FT_label = 30;

%% Load data

% load EnKF result
load EnKF_Lorenz96_data;

% load EnSRF result
load EnSRF_Lorenz96_Potter_data;

% load PF result
load PF_Lorenz96_data;

hfig = figure;
    set(hfig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% title(['Emsemble number = ', num2str(Nens)]);

NumPlot = 1;            % Number of subplots
for j = 1 : NumPlot
   
    subplot(NumPlot, 1, j);
    plot(tPredict, yPredict(:, Obs(j)),'r--',...
        tReference, yReference(:, Obs(j)),'k-.', ...
        tAnalysis_EnKF, yAnalysis_EnKF(:, Obs(j)),'b','MarkerSize', MS, 'LineWidth', 3);
%         tAnalysis_EnSRF_Potter, yAnalysis_EnSRF_Potter(:, Obs(j)),'g-',...
%         tAnalysis_PF, yAnalysis_PF(:, Obs(j)),'m-',...
%         tReference(ObsPoints), yReference(ObsPoints, Obs(j)), 'ko',...
        
    xlabel('Time steps', 'fontsize', FT_label, 'FontWeight','bold');
    ylabel('System State', 'fontsize', FT_label, 'FontWeight','bold');
    zoom on;
    set(gca,'FontSize',FS);
    if j==1
        h = legend('forecast', 'reference','EnKF-analysis');
        set(h,'FontSize',LF);
        legend boxoff;
    end
   
   
end

-- 
Haiyan Cheng
Associate Professor
Computer Science Department
Willamette University
900 State St.
Salem, Oregon, 97301
