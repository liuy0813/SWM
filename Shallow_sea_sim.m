function waterwave () %Drop_height, time

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

% Set up Reference Matrix
% fprintf('Loading reference matrix...\n')
% %ref_mat = zeros(100, 500,  64, 64, 3);
% mat = importdata('REF_matrix.mat','-mat'); % default var name is ref_mat
%
% [ens,time,Nvar,~,vals] = size(mat);  % read in matrix
fprintf('Bulding Obs matrix for H...\n')

ObsValuesH = importdata('OBS_matrix_H.mat','-mat');
% ObsValuesU = importdata('OBS_matrix_U.mat','-mat');
% ObsValuesV = importdata('OBS_matrix_V.mat','-mat');
[time,Nvar,~] = size(ObsValuesH);
% for i = 1 : time
%     ObsValuesH(i,:,:) = squeeze(mean(mat(:,i,:,:,1)));
%     %ObsValuesU(i,:,:) = squeeze(mean(mat(:,i,:,:,2)));
%     %ObsValuesV(i,:,:) = squeeze(mean(mat(:,i,:,:,3)));
% end

m=1; % for now vet state variables to 1, will be 3 in the end

% Observation mapping operator, H is a fat short matrix
Obs = 1 : Nvar; % Observe all
H_obs = zeros(m, (Nvar+2)^2);
for i = 1 : m
    H_obs(i, Obs(i)) = 1;
end

%% define conditions of ensamble
Nens = 10;
time = 50;

g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;         % plot interval
%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval

%% Initial Drop
drop_dim = 21;
D = zeros(21,21,Nens);  % create empty array for different drops

for i = 1 : Nens
    a = 1.0;                  % min size
    b = 3.0;                  % max size
    height = (b-a).*randn(1,1) + a;   % initial drop size
    D(:,:,i) = droplet(height,drop_dim);     % simulate a water drop (size,???)
end
%% Make empty vector for RMS of ens

RMS_ens = zeros(time,1);

%% Init. graphics
[surfplot,top] = initgraphics(Nvar);


tic;

%% Create ensamble of zeros/ones to store models

H = ones(Nvar+2,Nvar+2,Nens);   U = zeros(Nvar+2,Nvar+2,Nens);  V  = zeros(Nvar+2,Nvar+2,Nens);
Hx  = zeros(Nvar+1,Nvar+1,Nens); Ux  = zeros(Nvar+1,Nvar+1,Nens); Vx  = zeros(Nvar+1,Nvar+1,Nens);
Hy  = zeros(Nvar+1,Nvar+1,Nens); Uy  = zeros(Nvar+1,Nvar+1,Nens); Vy  = zeros(Nvar+1,Nvar+1,Nens);


%% Run Shallow Water Model
for itime = 1 : time
    % Debugging code!!!
    if mod(itime, 5) == 0
        fprintf('Current run: %d \n',itime)
    elseif itime <= 5
        fprintf('Current run: %d \n',itime)
    end
    
    % Observations at the end of the time interval
    z = (squeeze(ObsValuesH(itime, :,:)))';
    
    % measurement error and covariance
    gama = 0 + 1.*randn(m, Nens);       % Gaussian observation perturbation, Generate values from a normal distribution with mean 0 and standard deviation 1.
    Obs_cov = (gama * gama') / (Nens - 1);
    Obs_ens = zeros(m, Nens);
    
    % perturbed measurement
    for i = 1 : m
        Obs_ens(i, :) = z(i) + gama(i, :);    %% Measurement Ensemble
    end
    
    
    % initialize water drops
    for k = 1 : Nens
        if itime == 1;
            w = size(D(:,:,k),1);
            i = 5 +(1:w);
            j = 5 +(1:w);
            H(i,j,k) = H(i,j,k) + 0.5*D(:,:,k);
        end
        
        % Reflective boundary conditions
        H(:,1,k) = H(:,2,k);
        U(:,1,k) = U(:,2,k);
        V(:,1,k) = -V(:,2,k);
        H(:,Nvar+2,k) = H(:,Nvar+1,k);
        U(:,Nvar+2,k) = U(:,Nvar+1,k);
        V(:,Nvar+2,k) = -V(:,Nvar+1,k);
        H(1,:,k) = H(2,:,k);
        U(1,:,k) = -U(2,:,k);
        V(1,:,k) = V(2,:,k);
        H(Nvar+2,:,k) = H(Nvar+1,:,k);
        U(Nvar+2,:,k) = -U(Nvar+1,:,k);
        V(Nvar+2,:,k) = V(Nvar+1,:,k);
        
        %% Take a half time step to estimate derivatives at middle time.
        
        % x direction
        i = 1:Nvar+1;
        j = 1:Nvar;
        
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
        i = 1:Nvar;
        j = 1:Nvar+1;
        
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
        
        i = 2:Nvar+1;
        j = 2:Nvar+1;
        
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
        

        %% Compute ensemble
        
        OneN(1 : Nens, 1 : Nens) = 1 / Nens;
        H = permute(H,[1 2 3]);
        H = reshape(H,(Nvar+2)^2,10);
        Hbar   = H * OneN; % can't use three dimensions
        Hprime = H - Hbar;
        
        %% Compute ensemble covariance matrix
        Ens_cov = (Hprime * Hprime') / (Nens - 1);
        
        M = H_obs * Ens_cov * H_obs' + Obs_cov;   % Analysis equation
        
        %% Compute M inverse
        %handling the singular values
        
        [Uni_mat_U, S, Uni_mat_V] = svd(M);   % single value decomposition
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
        
        mInverse = Uni_mat_V * S * Uni_mat_U';    % M inverse = V * inv(S) * U';
        
        H = H + Ens_cov * H_obs' * ( mInverse * (Obs_ens - H_obs * H) );
        
        % Reshape H back to 3 dimensions!!!
        H = reshape(H,Nvar+2,Nvar+2,10);
        
        % Update plot
        if mod(k,Nens) == 0
            C = abs(U(i,j,k)) + abs(V(i,j,k));  % Color shows momemtum
            set(surfplot,'zdata',H(i,j,k),'cdata',C);
            set(top,'string',sprintf('step = %d',itime))
            drawnow
        end
    end
    
end
disp('Run time.....');
toc;

%% Save result for plotting and post-analysis purpose
filename = 'EnKF_SWM.mat';
save (filename);
end
% ------------------------------------

function D = droplet ( height, width )

% DROPLET  2D Gaussian
% D = droplet(height,width)
%
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
c = (1:64)'/64;
cyan = [c*0 c c];
colormap(cyan)
top = title('Shallow Sea Sim Ensemble');

return
end

%  GOES IN MAIN METHOD
%         if(mod(nstep,sample) == 0)
%             %make avg matrix for state - all 3 variables for all x,y
%             Tstart = tic
%             ens_avg = get_ens_avg(H,U,V);
%             obs = get_obs(nstep,0); %update to change error
%             t1 = toc(Tstart)
%
%             %generate errors
%             ens_err = get_ens_err(H,U,V,ens_avg,nstep);
%             obs_err = get_obs_err(0);
%             t2 = toc(Tstart)
%
%             %calc kalman gain
%             k_gain = get_k_gain(ens_err,obs_err);
%             t3 = toc(Tstart)
%             %ofset ens_avg
%             analysis_change = get_ana_chng(ens_avg,obs,k_gain);
%             t4 = toc(Tstart)
%
%             %update state
%             update_ens(analysis_change,H,U,V);
%             t2 = toc(Tstart)
%
%         end

% function ens_sum = get_ens_sum(H,U,V)
% sum_mat = zeros(64,64,3);
% ensemble = size(H,3);
% x = 1:64;
% y=1:64;
% xR = 2:65;
% yR=2:65;
%         for i=1:ensemble
%             sum_mat(x,y,1) = sum_mat(x,y,1) + H(xR,yR,i);
%             sum_mat(x,y,2) = sum_mat(x,y,2) + U(xR,yR,i);
%             sum_mat(x,y,3) = sum_mat(x,y,3) + V(xR,yR,i);
%         end
%
% ens_sum = sum_mat;
% end
%
% function avg = get_ens_avg(H,U,V)
% ensemble = size(H,3);
% avg = get_ens_sum(H,U,V) ./ ensemble;
% end
%
% function obs = get_obs(time,error)
% ensemble_num = 100; %MAGIC NUMBER
% sum = zeros(64,64,3);
% file = importdata('REF_matrix.mat','-mat');
%
% for x=1:64
%     for y=1:64
%         for data=1:3
%             for run = 1:ensemble_num
%                 sum(x,y,data) = sum(x,y,data)+file(run,time,x,y,data);
%             end
%         end
%     end
% end
%
% avg = sum ./ ensemble_num;
% %do errror
% obs = avg;
% end
%
% function error_ens = get_ens_err(H,U,V,ens_sum,T)
% total_err = zeros(64,64,3);
% ensemble = size(H,3);
% for i=1:ensemble
%     for x=1:64
%         for y = 1:64
%             total_err(x,y,1) = total_err(x,y,1) + (H(x,y,i) - ens_sum(x,y,1))...
%                 .*(H(x,y,1) - ens_sum(x,y,1))^T;
%
%             total_err(x,y,2) = total_err(x,y,2) + (U(x,y,i) - ens_sum(x,y,2))...
%                 .*(H(x,y,i) - ens_sum(x,y,2))^T;
%
%             total_err(x,y,3) = total_err(x,y,3) + (V(x,y,i) - ens_sum(x,y,3))...
%                 .*(H(x,y,i) - ens_sum(x,y,3))^T;
%         end
%     end
% end
% error_ens = total_err/(ensemble-1);
% end
%
% function obs_err = get_obs_err(input)
% obs_err = eye(3); %we don't know what to do here
% end
%
% function k_gain = get_k_gain(ens_err, obs_err)
% H = eye(3,1);
% k_gain = zeros(64,64,3);
% size(ens_err)
% size(obs_err)
% obs_err
% size(H)
% for i = 1:64
%     for j = 1:64
%         k_gain(i,j,:) = ens_err(i,j,:) .* (H.*ens_err(i,j,:)+obs_err)^-1; %what is H^t
%     end
% end
% end
%
% function analysis_change = get_ana_chng(ens_avg,obs,k_gain)
%
% analysis_change =k_gain(obs - ens_avg);
% end
%
% function update_ens(ana_chng,H,U,V)
% ensemble = size(H,3);
% for x=1:64
%     for y = 1:64;
%         for i=1:ensemble
%             H(x,y,i) = H(x,y,i) + ana_chng(x,y,1);
%             U(x,y,i) = U(x,y,i) + ana_chng(x,y,2);
%             V(x,y,i) = V(x,y,i) + ana_chng(x,y,3);
%         end
%     end
% end
% end

