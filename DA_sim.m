function DA_sim(Nens,time, deff_da_freq) %Drop_height, time
close all
clc

%% WATER WAVE
% 2D Shallow Water Model
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.
% Random water drops initiate gravity waves.
% Surface plot displays height colored by momentum.
% Plot title shows t = simulated time and tv = a measure of total variation.
% An exact solution to the conservation law would have constant tv.
% Lax-Wendroff produces nonphysical oscillations and increasing tv.
%
%  Author:
%
%    Cleve Moler
%
% Edited by:
%
% Carl Rodriguez, Tyler Welch
%(Willamette University, NSF research)
%
%  Reference:
%
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www


%% Deffault arguments
if nargin == 0
    Nens = 20;
    time = 500;
    deff_da_freq = 50;
end
fprintf('Running on %d ensemble size, for %d time\n',Nens,time)

%%
fprintf('Bulding Obs matrix for H...\n')
%% define conditions of ensamble


ObsValuesH = importdata('Data/OBS_matrix_H.mat','-mat');
% ObsValuesU = importdata('Data/OBS_matrix_U.mat','-mat');
% ObsValuesV = importdata('Data/OBS_matrix_V.mat','-mat');
[~,xDim,yDim] = size(ObsValuesH);

%% define model enviornment

g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;         % plot interval

%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval

%% Define DA constants
num_elems =xDim*yDim; % number of elements in grid
Obs = 1 : num_elems; % Observe all, can change step size to alter frequency of observations
H_Map = zeros(num_elems);
var_array = zeros(time,1);

for i = 1 : num_elems
    H_Map(i, Obs(i)) = 1;
end
%% Initial Drop
drop_dim = 21;

D = zeros(21,21,Nens);  % create empty array for different drops

for i = 1 : Nens
    center = 1.5;
    std_dev = .1;% max size
    height = center - 1  + std_dev*randn(1,1); % highest point of droplet
    %the initial drop is added onto water of height 1
    %so 1 is subtracted from center
    
    D(:,:,i) = droplet(height,drop_dim); %create gaussian droplet
end

%% Make empty vector for RMSE of EnKF vs. REF



RMSE = zeros(time,1);
RMS_H = zeros(time,xDim,yDim);
temp = 0;
temp2 = zeros(Nens,xDim,yDim);
% test_H_mean = zeros(time,1);
% test_H = zeros(time,1);


%% variables for pdfs

% Nbins = 10;
DA_num = 0;
curr_da_freq = deff_da_freq;
DA_deff_size = floor(time/deff_da_freq);
DA_curr_size = DA_deff_size;
DA_loc = zeros(1,DA_deff_size);
avg_div = zeros(DA_deff_size,1);

markers = 1:deff_da_freq:time;


pdf_coords = zeros(DA_deff_size, xDim, yDim, Nens);
distance_R = zeros(DA_deff_size,1);
distance_L = zeros(DA_deff_size,1);
hist_size = Nens;
waterfall_data_pre = zeros(DA_deff_size,hist_size);
waterfall_data_post = zeros(DA_deff_size,hist_size);
% hist_prior = zeros(time,xDim,yDim,Nens);
% hist_post= zeros(time,xDim,yDim,Nens);

pdfs_prior = zeros(DA_deff_size,xDim,yDim,hist_size);
pdfs_post = zeros(DA_deff_size,xDim,yDim,hist_size);


%% Init. graphics
% [surfplot,top] = initgraphics(xDim);

%% Init. timer

tic;

%% Create ensamble of zeros/ones to store models

H = ones(xDim+2,yDim+2,Nens);   U = zeros(xDim+2,yDim+2,Nens);  V  = zeros(xDim+2,xDim+2,Nens);
Hx  = zeros(xDim+1,yDim+1,Nens); Ux  = zeros(xDim+1,yDim+1,Nens); Vx  = zeros(xDim+1,xDim+1,Nens);
Hy  = zeros(xDim+1,yDim+1,Nens); Uy  = zeros(xDim+1,yDim+1,Nens); Vy  = zeros(xDim+1,xDim+1,Nens);


%% Run Shallow Water Model
for itime = 1 : time
    %% Output first 5 loops
    if mod(itime, 5) == 0
        fprintf('Current run: %d \n',itime)
        % Increase loop message frequency
    elseif itime <= 10
        fprintf('Current run: %d \n',itime)
    end
    
    
    
    
    max_vals = zeros(Nens,1);
    % Nested loop for Ensembles
    for k = 1 : Nens
        
        %% initialize water drop
        if itime == 1;
            w = size(D(:,:,k),1);
            i = 5 +(1:w);
            j = 5 +(1:w);
            H(i,j,k) = H(i,j,k) + D(:,:,k);
            max_vals(k) = max(max(H(:,:,k)));
        end
        
        %% Reflective boundary conditions
        H(:,1,k) = H(:,2,k);
        U(:,1,k) = U(:,2,k);
        V(:,1,k) = -V(:,2,k);
        H(:,yDim+2,k) = H(:,yDim+1,k);
        U(:,yDim+2,k) = U(:,yDim+1,k);
        V(:,yDim+2,k) = -V(:,yDim+1,k);
        H(1,:,k) = H(2,:,k);
        U(1,:,k) = -U(2,:,k);
        V(1,:,k) = V(2,:,k);
        H(xDim+2,:,k) = H(xDim+1,:,k);
        U(xDim+2,:,k) = -U(xDim+1,:,k);
        V(xDim+2,:,k) = V(xDim+1,:,k);
        
        %% Take a half time step to estimate derivatives at middle time.
        
        % x direction
        i = 1:xDim+1;
        j = 1:yDim;
        
        % height
        Hx(i,j,k) = (H(i+1,j+1,k)+H(i,j+1,k))/2 - dt/(2*dx)*(U(i+1,j+1,k)-U(i,j+1,k));
        
        % x momentum
        Ux(i,j,k) = (U(i+1,j+1,k)+U(i,j+1,k))/2 -  ...
            dt/(2*dx)*((U(i+1,j+1,k).^2./H(i+1,j+1,k) + g/2*H(i+1,j+1,k).^2) - ...
            (U(i,j+1,k).^2./H(i,j+1,k) + g/2*H(i,j+1,k).^2));
        
        %             % y momentum
        Vx(i,j,k) = (V(i+1,j+1,k)+V(i,j+1,k))/2 - ...
            dt/(2*dx)*((U(i+1,j+1,k).*V(i+1,j+1,k)./H(i+1,j+1,k)) - ...
            (U(i,j+1,k).*V(i,j+1,k)./H(i,j+1,k)));
        
        % y direction
        i = 1:xDim;
        j = 1:yDim+1;
        
        % height
        Hy(i,j,k) = (H(i+1,j+1,k)+H(i+1,j,k))/2 - dt/(2*dy)*(V(i+1,j+1,k)-V(i+1,j,k));
        
        % x momentum
        Uy(i,j,k) = (U(i+1,j+1,k)+U(i+1,j,k))/2 - ...
            dt/(2*dy)*((V(i+1,j+1,k).*U(i+1,j+1,k)./H(i+1,j+1,k)) - ...
            (V(i+1,j,k).*U(i+1,j,k)./H(i+1,j,k)));
        %             % y momentum
        Vy(i,j,k) = (V(i+1,j+1,k)+V(i+1,j,k))/2 - ...
            dt/(2*dy)*((V(i+1,j+1,k).^2./H(i+1,j+1,k) + g/2*H(i+1,j+1,k).^2) - ...
            (V(i+1,j,k).^2./H(i+1,j,k) + g/2*H(i+1,j,k).^2));
        
        %% Now take a full step that uses derivatives at middle point.
        
        i = 2:xDim+1;
        j = 2:yDim+1;
        
        % height
        H(i,j,k) = H(i,j,k) - (dt/dx)*(Ux(i,j-1,k)-Ux(i-1,j-1,k)) - ...
            (dt/dy)*(Vy(i-1,j,k)-Vy(i-1,j-1,k));
        %  x momentum
        U(i,j,k) = U(i,j,k) - (dt/dx)*((Ux(i,j-1,k).^2./Hx(i,j-1,k) + g/2*Hx(i,j-1,k).^2) - ...
            (Ux(i-1,j-1,k).^2./Hx(i-1,j-1,k) + g/2*Hx(i-1,j-1,k).^2)) ...
            - (dt/dy)*((Vy(i-1,j,k).*Uy(i-1,j,k)./Hy(i-1,j,k)) - ...
            (Vy(i-1,j-1,k).*Uy(i-1,j-1,k)./Hy(i-1,j-1,k)));
        %  y momentum
        V(i,j,k) = V(i,j,k) - (dt/dx)*((Ux(i,j-1,k).*Vx(i,j-1,k)./Hx(i,j-1,k)) - ...
            (Ux(i-1,j-1,k).*Vx(i-1,j-1,k)./Hx(i-1,j-1,k))) ...
            - (dt/dy)*((Vy(i-1,j,k).^2./Hy(i-1,j,k) + g/2*Hy(i-1,j,k).^2) - ...
            (Vy(i-1,j-1,k).^2./Hy(i-1,j-1,k) + g/2*Hy(i-1,j-1,k).^2));
    end
    
    
    
    
    %%Determine if DA is performed
    do_DA = false;
    
    if DA_num > 0 && itime ==  DA_loc(1,DA_num) + curr_da_freq
        do_DA = true;
    elseif DA_num == 0 && itime == curr_da_freq
        do_DA = true;
    end
    
    if do_DA
        %%Perform DA
        DA_num = DA_num + 1;
        if DA_num > DA_curr_size
            DA_curr_size = DA_curr_size + 1;
            %maybe reallocate
        end
        DA_loc(1,DA_num) = itime;
        
        tic;
        fprintf('Starting DA at time: %d \n', itime)
        % Observations at the end of the time interval
        
        
        z = squeeze(ObsValuesH(itime, :,:))'; %check itime is properly done
        zsq = squeeze(reshape(z,xDim*yDim,1));
        
        % measurement error and covariance
        error = 0.05;
        gama = repmat(error, num_elems, Nens);   % Gaussian observation perturbation, Generate values from a normal distribution with mean 0 and standard deviation 1.
        
        
        Obs_cov = (gama * gama') / (Nens - 1);
        Obs_ens = zeros(num_elems, Nens);
        
        % perturbed measurement
        for ens = 1 : Nens
            Obs_ens(:,ens) = zsq + error*randn(num_elems,1);    %% Measurement Ensemble
        end
        
        % Reshape data into one column
        
        OneN(1 : Nens, 1 : Nens) = 1 / Nens;
        Hpre = H;
        Hperm = permute(H,[1 2 3]);
        Hresh = reshape(Hperm(2:xDim+1,2:yDim+1,:),xDim*yDim,Nens);
        
        Hbar   = Hresh * OneN; % can't use three dimensions
        Hprime = Hresh - Hbar;
        
        % Compute ensemble covariance matrix
        Ens_cov = (Hprime * Hprime') / (Nens - 1);
        
        M = Ens_cov + Obs_cov;   % Analysis equation
        
        %Do svd calc for inverse - takes most time of all assimilation
        fprintf('Beginning svd calc \n')
        tic;
        
        [Uni_mat_U, S, Uni_mat_V] = svd(M);   % single value decomposition
        
        
        Xi = diag(S);% singular values
        
        nos = length(Xi);
        
        for jj = 1 : nos
            if (Xi(jj) < Xi(1) / (10^6))
                Xi(jj) = 0;
            else
                Xi(jj) = 1 / Xi(jj);
            end
        end
        
        S = diag(Xi);
        
        
        fprintf('SVD time: %3.1f seconds. Starting mInverse \n',toc)
        
        tic;
        mInverse = Uni_mat_V * S * Uni_mat_U';    % M inverse = V * inv(S) * U';
        
        fprintf('Done with inverse, creating analysis\n')
        %% Data Assimilate
        %         diff = Obs_ens - H_Map * Hresh;
        
        % Inflation code
        %         inf_factor = 1.5;
        %         Ens_cov = Ens_cov * inf_factor;
        
        Hresh = Hresh + Ens_cov * H_Map' * ( mInverse * (Obs_ens - H_Map * Hresh) );
        %         mInverse;
        %         meaninverse = mean(mean(mInverse))
        %         rangeinverse = mean(range(mInverse))
        %         Hresh = Hresh + 0.2*(Obs_ens - H_Map * Hresh);
        
        fprintf('Time for inverse+analysis: %d reshaping\n',toc)
        
        
        % Reshape H back to 3 dimensions!!!
        H = reshape(Hresh,xDim,yDim,Nens);
        Obs_ens_resh = reshape(Obs_ens,xDim,yDim,Nens);
        
        
        B = zeros(xDim+2,yDim+2,Nens);
        B(2:xDim+1,2:yDim+1,:)= H;
        H = B;
        x = 16;
        y = 16;
        fprintf(' pre: %2.5f obs_ens: %2.5f  post: %2.5f  \n',...
            Hpre(x,y,:),Obs_ens_resh(x,y,:),H(x,y,:))
        fprintf('obs values: %2.5f \n',ObsValuesH(itime,x-1,y-1))
        fprintf('Done with DA. \n')
        
        %% Create normal distribution of pre and post
        sum = 0;
        fprintf('Starting to create normal histograms')
        for q = 1:xDim
            for p = 1:xDim
                %                 hist_post(itime,q,p,:) = hist(H(q+1,p+1,:)); % hist on third dim, i.e. ensembles
                minV = min(min(H(q+1,p+1,:),Hpre(q+1,p+1,:)));
                maxV = max(max(H(q+1,p+1,:),Hpre(q+1,p+1,:)));
                
                x_vals = linspace(minV,maxV,hist_size);
                
                % Store pdf x coordinates to plot on singple graph
                
                pd1 = fitdist(squeeze(Hpre(q+1,p+1,:)), 'Normal');
                
                pdfs_prior(DA_num,q,p,:) = pdf(pd1,x_vals);
                
                pd2 = fitdist(squeeze(H(q+1,p+1,:)), 'Normal');
                pdfs_post(DA_num,q,p,:) = pdf(pd2,x_vals);
                
                if maxV ~= minV
                    pdf_coords(DA_num,q,p,:) = minV:(maxV-minV)/(Nens-1):maxV;
                end
                x1 = squeeze((pdfs_post(DA_num,x,y,:)))';
                x2 = squeeze((pdfs_prior(DA_num,x,y,:)))';
                if any(x1) > 0 && any(x2) > 0
                    div = KLDiv(x2,x1);
                    sum = sum + div;
                end
                
                if p == 15 && q == 15
                    fprintf('Range pre: %d, range post %d',....
                        range(Hpre(q+1,q+1,:)),range(H(q+1,q+1,:)))
                    squeeze(Hpre(q+1,p+1,:))
                    squeeze(H(q+1,p+1,:))
                    x_vals_location = x_vals;
                    size(minV:(maxV-minV)/(Nens-1):maxV)
                    
                end
                
            end
        end
        
        %% Plot pre and post distributions
        
        figure(1)
        waterfall_data_pre(DA_num,:) = squeeze(pdfs_prior(DA_num,16,16,:));
        waterfall_data_post(DA_num,:) = squeeze(pdfs_post(DA_num,16,16,:));
        
        x_vals
        squeeze(pdfs_prior(DA_num,16,16,:))
        squeeze(pdfs_post(DA_num,16,16,:))
        
        plot( x_vals_location,squeeze(pdfs_prior(DA_num,16,16,:)),'b', x_vals_location,squeeze(pdfs_post(DA_num,16,16,:)),'r', 'LineWidth', 2)
        legend('PDF prior DA','PDF post DA')
        title('Probability Density function before and after DA', 'fontsize', 20, 'fontweight', ...
            'bold');
        ylabel('probability', 'fontsize', 15, 'fontweight', 'bold');
        xlabel('height at point 16,16', 'fontsize', 15, 'fontweight', 'bold');
        name=['Data/fig',num2str(DA_num),'.png'];
        %saveas(gca,name);
        %         figure('units','normalized','outerposition',[0 0 1 1])
        
        avg_div(DA_num) = sum/(xDim*yDim);
        avg_div(DA_num)
        if avg_div(DA_num) > 0.5
            curr_da_freq = ceil(curr_da_freq/2)
        elseif avg_div(DA_num) < 0.4
            curr_da_freq = ceil(curr_da_freq*2)
        end
        
        
        
        % end of DA
    end
    %% Calc Error
    for q = 1:xDim
        for p = 1:xDim
            for s = 1 : Nens
                temp = temp +  (H(q+1,p+1,s)-ObsValuesH(itime,q,p))^2;
                temp2(s,q,p) = H(q+1,p+1,s);
            end
            RMS_H(itime,q,p) = sqrt(temp/Nens);
            temp = 0;
        end
    end
    
    var_array(itime) = mean(mean(var(temp2)));
    
    %% Update plot
    i = 2:xDim+1;
    j = 2:yDim+1;
    %     test_H_mean(itime)  = mean(H(16,16,:));
    %     test_H(itime)  = H(16,16,Nens);
    %
    %     C = abs(U(i,j,Nens)) + abs(V(i,j,Nens));  % Color shows momemtum
    %     set(surfplot,'zdata',H(i,j,Nens),'cdata',C);
    %     set(top,'string',sprintf('step = %d',itime))
    drawnow
    
    
    %% Check distribution
    
    if(false && itime == 1)%true to check variation of initial droplet
        check_std_dev(std_dev,center,max_vals)
    end
end
disp('Run time.....');
toc;
avg_div


for i = 1 : time
    RMSE(i) =  mean(mean(RMS_H(i,:,:)));
    % RMSE(i) =  RMS_H(i,16,16);
end

count = 1;
for i = 1:DA_num
    x1 = squeeze((pdfs_post(i,16,16,:)))';
    x2 = squeeze((pdfs_prior(i,16,16,:)))';
    
    distance_R(count) = KLDiv(x1,x2);
    distance_L(count) = KLDiv(x2,x1);
    count = count + 1;
end


%% Plot results and error

%add error bars for variance!!! var of ensemble
% figure(2)
% size(RMSE)
% plot(xtime,RMSE,'g','LineWidth',1)
% title('RMS Error of Reference Model', 'fontsize', 20, 'fontweight', ...
%     'bold');
% xlabel('time', 'fontsize', 15, 'fontweight', 'bold');
% ylabel('RMSE', 'fontsize', 15, 'fontweight', 'bold');
% 
% Obs_point = ObsValuesH(xtime,16,16);
% xtime = 1:time;


y_coords = zeros(DA_num,Nens);
for i = 1 : DA_num
    y_coords(i,:) = i;
end

% C1 = colormap(hsv);
% C2 = colormap(bone);

size(pdf_coords)
size(y_coords)
size(squeeze(pdfs_prior(:,16,16,:)))


figure(4)
for q = 1:xDim
    for p = 1:xDim
        if max(squeeze(pdfs_prior(:,q,p,:))) < 100
            waterfall(squeeze(pdf_coords(:,q,p,:)),y_coords,squeeze(pdfs_prior(:,q,p,:)))
            hold on
        end
    end
end
% waterfall(pdf_coords,y_coords,squeeze(pdfs_post(:,16,16,:)))
title('Pdfs Prior and Post DA', 'fontsize', 20, 'fontweight', ...
    'bold')
xlabel('height', 'fontsize', 15, 'fontweight', 'bold');
ylabel('time', 'fontsize', 15, 'fontweight', 'bold');
zlabel('probability', 'fontsize', 15, 'fontweight', 'bold');
hold off

% figure(5)
% waterfall(pdf_coords,y_coords,squeeze(pdfs_post(:,16,16,:)))
% legend('Height','Observation', 'Probability')
% title('Pdfs post')

figure(6)
DA_loc
plot(DA_loc,distance_R,'b-*',DA_loc,distance_L,'r-+','LineWidth',2,'LineSmoothing','on')
title('Distance between hists prior and post DA', 'fontsize', 20, 'fontweight', ...
    'bold');
legend('Distance Right', 'Distance Left');
xlabel('time', 'fontsize', 15, 'fontweight', 'bold');
ylabel('distance/divergence', 'fontsize', 15, 'fontweight', 'bold');

%% Save result for plotting and post-analysis purpose
filename = 'Data/EnKF_SWM.mat';
save (filename);
save('Data/EnKF_error.mat','RMSE');
save('Data/drops.mat', 'D');
save('Data/var_EnKF.mat','var_array');
end
% ------------------------------------

function D = droplet ( height, width )

% DROPLET  2D Gaussian
% D = droplet(height,width)
[ x, y ] = ndgrid ( -1:(2/(width-1)):1 );

D = height * exp ( -5 * ( x.^2 + y.^2 ) );

return
end

% ------------------------------------
function [surfplot,top] = initgraphics(n)

% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

clf
shg
set(gcf,'numbertitle','off','name','Shallow_water')
x = (0:n-1)/(n-1);
surfplot = surf(x,x,ones(n,n),zeros(n,n));
grid off
axis([0 1 0 1 -1 3])
caxis([-1 1])
shading faceted
c = (1:64)'/64; %not sure if this should be x or y dim
cyan = [c*0 c c];
colormap(winter)
top = title('Shallow Sea Sim Ensemble');

return
end

%----------------------------------------
function KL = getKLD(Hpre,H)
%bins = numbe of ensembles / 3
xDim = size(H,1)-2;
yDim = size(H,2)-2;
KLs = zeros(xDim,yDim);
% for x = 2 : 3
%     for y = 2 : 3
Hpre(16,16,:)
H(16,16,:)
a = makeHist(Hpre(16,16,:),H(16,16,:));
b = makeHist(H(16,16,:),Hpre(16,16,:));
a
b
KLDiv(a,b)
KLs(16-1,16-1) = squeeze(KLDiv(b,a));
%     end
% end
KL = mean(mean(KLs));
end

function histo = makeHist(a,b)
%% create histogram of a based of data from a and b
a = squeeze(a);
b = squeeze(b);
assert(size(a,1) == size(b,1),'Input size does not match in makeHist');


minV = min(min(a,b))
maxV = max(max(a,b))

Nens = size(a,1);
bin_num = floor(Nens/3); % if we have fewer bins the data will be quite spread
assert( bin_num >= 3,'Below minimum number in ensemble to make histogram');
%we may want to alter this value
if minV == maxV
    fprintf('Making histogram with only one point')
    histo = ones(1,bin_num);
    return
end


bin_size = (maxV-minV)/bin_num;

assert(bin_size > 0);

histA = zeros(1,bin_num);

for i = 1:Nens
    diff = a(i) - minV;
    loc = ceil(diff/bin_size); %not using ceil for elem = minV
    if loc == 0
        assert(a(i) == minV,'Got a loc = 0 for something not minV');
        loc = 1;
    end
    histA(1,loc) = histA(1,loc) + 1;
end


histo = histA;

end

%-------------------------------------
function new_freq = get_da_freq(pdf_post, pdf_prior,DA_num, curr_da_freq)
pdf_post(DA_num,:,:,:)

end
%-------------------------------------
function check_std_dev(std_dev, center, max_vals)
fprintf('Max is %d and min is %d \n',max(max_vals),min(max_vals))
one_dev = 0;
two_dev = 0;
three_dev = 0;
for itr = 1:size(max_vals,1)
    if std_dev > abs(center - max_vals(itr))
        one_dev = one_dev +1;
    elseif std_dev > .5 * abs(center - max_vals(itr))
        two_dev = two_dev + 1;
    elseif std_dev > .25 *abs(center - max_vals(itr))
        three_dev = three_dev + 1;
    end
end
two_dev = two_dev+one_dev;
three_dev = two_dev+three_dev;
fprintf('One %d , two %d , three %d \n', one_dev, two_dev, three_dev);
end