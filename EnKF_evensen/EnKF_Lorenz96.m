%% EnKF applied to Lorenz-96 system

% function [R1, R2, RMS1, RMS2] = EnKF_Lorenz96(Nens)

%% Initialize model parameters
initialization;

disp('EnKF ensemble.....');
disp(['After ', num2str(spinup), ' spinup with ', ...
    num2str(Nens), ' ensemble members.']);

tic;

%% Apply EnKF with perturbed observations for each ensemble member as described in Evensen paper.

% Integration from one obs. to another obs. point

for itime = 2 : length(ObsPoints)
    % itime
    % Time interval
    % Tinterval = [ObsTimes(itime-1),ObsTimes(itime)];
    Tinterval1 = ObsTimes(itime - 1) : 0.01 : ObsTimes(itime);
    Tgrid = linspace(ObsTimes(itime - 1), ObsTimes(itime), 10);
    
    % Observations at the end of the time interval
    z = ObsValues(itime, :)';
    
    % measurement error and covariance
    gama = randn(m, Nens);       % Gaussian observation perturbation 
    Obs_cov = (gama * gama') / (Nens - 1);
    
    Obs_ens = zeros(m, Nens);
    % perturbed measurement
    for i = 1 : m
        Obs_ens(i, :) = z(i) + gama(i, :);    %% Measurement Ensemble
    end
    
    %%  Integration of ensemble members
    
    for i = 1 : Nens
        start = A(:, i);     % Take one ensemble member
        [ti, yi] = ode45('lorenz96', Tinterval1, start);
        A(:, i) = yi(length(ti), :);
        tt = ti(length(ti));
        
        % Interpolation of the intermediate points
        for j = 1 : Nvar
            Y(:, j, i) = spline(ti, yi(:, j), Tgrid);
        end
    end
    
    %% Compute ensemble 
    
    OneN(1 : Nens, 1 : Nens) = 1 / Nens;
    Abar   = A * OneN;
    Aprime = A - Abar;
    
    %% Compute ensemble covariance matrix
    Ens_cov = (Aprime * Aprime') / (Nens - 1);
    
    % Perform Analysis
    % A=A+Pe*H'*( (H*Pe*H'+Re)\(D-H*A) );
    % M = H * Pe * H'+ 0.01 * H * Pe * H';
    
    M = H * Ens_cov * H' + Obs_cov;   % HPH'+R term
    
    %% Compute M inverse
    %handling the singular values
    
    [U, S, V] = svd(M);   % M=USV'
    Xi = diag(S);    % singular values
    nos = length(Xi);
    
    for jj = 1 : nos
        if (Xi(jj) < Xi(1) / (10^6))
            Xi(jj) = 0;
        else
            Xi(jj) = 1 / Xi(jj);
        end       
    end
    
    S = diag(Xi);
    
    mInverse = V * S * U';    % M inverse = V * inv(S) * U';
     
    A = A + Ens_cov * H' * ( mInverse * (Obs_ens - H * A) );
    
    me = mean(Y, 3);
    
    tAnalysis = [tAnalysis; Tgrid'] ;
    yAnalysis = [yAnalysis; me];   
    yObservation = [yObservation; yi(length(ti), :)];
end

disp('Run time.....');
toc; 

%% Save result for plotting and post-analysis purpose

tAnalysis_EnKF = tAnalysis;
yAnalysis_EnKF = yAnalysis;

filename = 'EnKF_Lorenz96_data.mat';
save (filename);

%% Plot results

Plot;


%% 


% figure
% 
% % title(['Emsemble number = ', num2str(Nens)]);
% NumPlot = 3;
% for j = 1 : NumPlot
%     subplot(NumPlot, 1, j);
%     plot(tp_EnKF, yp_EnKF(:, Obs(j)),'r--',...
%         tref, yref(:, Obs(j)),'k-.', ...
%         ta_EnKF, ya_EnKF(:, Obs(j)),'b',...
%         tref(ObsPoints), yref(ObsPoints, Obs(j)), 'ko',...
%         'MarkerSize', MS, 'LineWidth', 3);
%     set(gca,'FontSize',FS);
%     if j==1
%         h = legend('forecast', 'reference','analysis');
%         set(h,'FontSize',LF);
%         % legend boxoff;
%     end
% end
% 
% % % plot first 3 Un-observed species
% % figure
% % for j = 1 : 2
% %     subplot(2, 1, j);
% %     plot(tp_EnKF, yp_EnKF(:, UnObs(j)), 'r--',...
% %         tref, yref(:, UnObs(j)), 'k-.', ...
% %         ta_EnKF, ya_EnKF(:, UnObs(j)), 'b',...
% %         tref(ObsPoints), yref(ObsPoints, UnObs(j)), 'ko',...
% %         'MarkerSize', MS, 'LineWidth', 3);
% %     set(gca,'FontSize',FS);
% %     if j==1
% %         h = legend('forecast', 'reference','analyzed');
% %         set(h,'FontSize',LF);
% %         % legend boxoff;
% %     end
% % end
% 
% % computing solution accuracy at the observation points.
% % ref--true solution
% % bck--predicted solution
% % ana--analyzed solution
% 
% tr = yref(ObsPoints, :) * H';
% pred = yp_EnKF(ObsPoints, :) * H';
% ana = yobs * H';
% 
% ref = yref(ObsPoints, :) * H';
% bck = yp_EnKF(ObsPoints, :) * H';
% ana = yobs * H';
% 
% [points, l2] = size(tr);
% 
% % rms1=zeros(points);
% 
% for i =1: points
%     rms_KF_BckRef(i) = norm(bck(i, :) - ref(i, :)) / sqrt(Nvar);
%     rms_KF_AnaRef(i) = norm(ana(i, :) - ref(i, :)) / sqrt(Nvar);
% end
% 
% for i =1: points
%     re_KF_BckRef(i) = norm(bck(i, :) - ref(i, :)) ./ norm(ref(i, :));
%     re_KF_AnaRef(i) = norm(ana(i, :) - ref(i, :)) ./ norm(ref(i, :));
% end
% 
% 
% [l1, l2] = size(tr);  % obs * Nvar
% NN = l1 * l2;
% p1 = reshape(tr, NN, 1);
% p2 = reshape(pred, NN, 1);
% p3 = reshape(ana, NN, 1);
% 
% % true solution vs. predicted solution
% 
% s1 = (NN * sum(p1 .* p2) - sum(p1) * sum(p2))^2;
% s2 = NN * (sum(p1 .^ 2)) - (sum(p1))^2;
% s3 = NN * (sum(p2 .^ 2)) - (sum(p2))^2;
% 
% R1 = s1 / (s2 * s3);
% 
% RMS1 = sqrt(1 / NN * sum((p1 - p2) .^2));
% 
% % true solution vs. analyzed solution.
% 
% t1 = (NN * sum(p1 .* p3) - sum(p1) * sum(p3))^2;
% t2 = NN * (sum(p1 .^ 2)) - (sum(p1))^2;
% t3 = NN * (sum(p3 .^ 2)) - (sum(p3))^2;
% 
% R2 = t1 / (t2 * t3);
% 
% RMS2 = sqrt(1 / NN * sum((p1 - p3) .^ 2));
% 





% end
