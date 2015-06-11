function waterwave (Nens,time, obs_freq) %Drop_height, time
%
%  clear all
% close all
clc
% WATER WAVE
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
nargin

if (isempty(Nens) || isempty(time))
    Nens = 5;
    time = 100;
end

fprintf('Bulding Obs matrix for H...\n')
%% define conditions of ensamble


ObsValuesH = importdata('OBS_matrix_H.mat','-mat');
% ObsValuesU = importdata('OBS_matrix_U.mat','-mat');
% ObsValuesV = importdata('OBS_matrix_V.mat','-mat');
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
H_Map = zeros(num_elems, num_elems);

for i = 1 : num_elems
    H_Map(i, Obs(i)) = 1;
end
%% Initial Drop
drop_dim = 21;

D = zeros(21,21,Nens);  % create empty array for different drops

for i = 1 : Nens
    center = 1.5;
    std_dev = .3;% max size
    height = center + 1.75*std_dev*randn(1,1);   % initial value for droplet
    %note: 1.75 is a constant that makes final drop size close to normally
    %distributed - determined by testing different constants on 10,000 ens
    D(:,:,i) = droplet(height,drop_dim);     % simulate a water drop (size,???)
end

%% Make empty vector for RMSE of EnKF vs. REF

RMSE = zeros(time,1);
RMS_H = zeros(time,xDim,yDim);
temp = 0;
test_H_mean = zeros(time,1);
test_H = zeros(time,1);

%% Init. graphics
[surfplot,top] = initgraphics(xDim);


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
            H(i,j,k) = H(i,j,k) + 0.5*D(:,:,k);
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
    
    %% Perform DA every obs_freq runs
    if itime == 50 %mod(itime,obs_freq) == 0
        tic;
        fprintf('Starting DA at time: %d \n', itime)
        % Observations at the end of the time interval

        z = squeeze(ObsValuesH(floor(itime/5), :,:))';
        z = squeeze(reshape(z,xDim*yDim,1));
        
        % measurement error and covariance
        gama = zeros(num_elems) + .01*randn(num_elems);   % Gaussian observation perturbation, Generate values from a normal distribution with mean 0 and standard deviation 1.
        gama = squeeze(gama);
        Obs_cov = (gama * gama') / (Nens - 1);
        Obs_ens = zeros(num_elems, Nens);
        
        
        % perturbed measurement
        for i = 1 : Nens
            Obs_ens(:,i) = zsq + gama(1);    %% Measurement Ensemble
        end
        
        
        % Compute ensemble
        
        OneN(1 : Nens, 1 : Nens) = 1 / Nens;
        Hpre = H;
        Hperm = permute(H,[1 2 3]);
        Hresh = reshape(Hperm(2:xDim+1,2:yDim+1,:),xDim*yDim,Nens);
        
        Hbar   = Hresh * OneN; % can't use three dimensions
        Hprime = Hresh - Hbar;
        
        % Compute ensemble covariance matrix
        Ens_cov = (Hprime * Hprime') / (Nens - 1);
        
        M = Ens_cov + Obs_cov;   % Analysis equation
        
        
        
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

        
        fprintf('SVD time: %d calc starting mInverse \n',toc)
        tic;
        mInverse = Uni_mat_V * S * Uni_mat_U';    % M inverse = V * inv(S) * U';
        
        fprintf('Done with inverse, creating analysis\n')
        % Data Assimilate
        Hresh = Hresh + Ens_cov * H_Map' * ( mInverse * (Obs_ens - H_Map * Hresh) );
        %H = H + (Obs_ens - H_Map * H);
        fprintf('Time for inverse+analysis: %d reshaping\n',toc)
        tic;
        
        % Reshape H back to 3 dimensions!!!
        H = reshape(Hresh,xDim,yDim,Nens);
        
        B = zeros(xDim+2,yDim+2,Nens);
        B(2:xDim+1,2:yDim+1,:)= H;
        H = B;
        asdf = H - Hpre;
        for x = 1:xDim+2
            for y = 1:yDim+2
                if asdf > 1
                    fprintf('x: %d, y: %d',asdf(x,y))
                end
            end
        end
        fprintf('Done with DA. Time to end: %d \n',toc)
        
    end
    %% Calc Error
    for q = 1:xDim
        for p = 1:xDim
            for s = 1 : Nens
                %                     if abs((H(q+1,p+1,s)-ObsValuesH(itime,q,p))) > 0.5
                %                         (H(q+1,p+1,s)-ObsValuesH(itime,q,p))
                %                     end
                temp = temp +  (H(q+1,p+1,s)-ObsValuesH(floor(itime/5),q,p))^2;
            end
            RMS_H(itime,q,p) = sqrt(temp/Nens);
            temp = 0;
        end
    end
        %% Update plot
    
    test_H_mean(itime)  = mean(H(16,16,:));
    test_H(itime)  = H(16,16,Nens);
    
    C = abs(U(i,j,Nens)) + abs(V(i,j,Nens));  % Color shows momemtum
    set(surfplot,'zdata',H(i,j,Nens),'cdata',C);
    set(top,'string',sprintf('step = %d',itime))
    drawnow
    
    %% Check distribution
    
    if(true & itime == 1)
        check_std_dev(std_dev,center,max_vals)
    end 
end
disp('Run time.....');
toc;

%% Save result for plotting and post-analysis purpose
filename = 'EnKF_SWM.mat';
save (filename);
xtime = 1:time;

for i = 1 : time
    RMSE(i) =  mean(mean(RMS_H(i,:,:)));
    % RMSE(i) =  RMS_H(i,16,16);
end

%file = 'EnKF_error.mat';
save('EnKF_error.mat','RMSE');
save('drops.mat', 'D');

%% Plot results and error
figure(1)
size(RMSE)
plot(xtime,RMSE,'g','LineWidth',1)
%axis([0 time -0.2 0.5])
title('mean of RMS Error of Reference Model')
xlabel('time')
ylabel('RMSE')

Obs_point = ObsValuesH(xtime,16,16);

ymarkers = 1:obs_freq:time;

figure(2)
plot(xtime,test_H_mean,'b',xtime,Obs_point,'r') %,'LineWidth',1 ymarkers,test_H_mean(ymarkers),'m'
legend('Mean Height','Observed Height')
legend
axis([0 time 0 5])
title('ENKF Compare ObsValuesH to H_mean')
xlabel('time')
ylabel('value')

figure(3)
plot(xtime,test_H,'b',xtime,Obs_point,'r') % ,'LineWidth',1
legend('Rand Ens Height','Observed Height')
axis([0 time 0 5])
title('ENKF Compare ObsValuesH to some ens of H')
xlabel('time')
ylabel('value')
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