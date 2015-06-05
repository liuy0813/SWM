function waterwave (Nens,time) %Drop_height, time

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

if (isempty(Nens) || isempty(time))
    Nens = 5;
    time = 100;
end

fprintf('Bulding Obs matrix for H...\n')
%% define conditions of ensamble


ObsValuesH = importdata('OBS_matrix_H.mat','-mat');
RMSE_EnKF = importdata('EnKF_error.mat','-mat');
D = importdata('drops.mat','-mat');
% ObsValuesU = importdata('OBS_matrix_U.mat','-mat');
% ObsValuesV = importdata('OBS_matrix_V.mat','-mat');
[~,Nvar,~] = size(ObsValuesH);

RMSE = zeros(time,1);
RMS_H = zeros(time,64,64);
temp = 0;

%% define model enviornment
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;         % plot interval
%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval

%% Initial Drop
% drop_dim = 21;
% D = zeros(21,21,Nens);  % create empty array for different drops
% 
% for i = 1 : Nens
%     a = 1.5;                  % min size
%     b = 2.5;                  % max size
%     height = (b-a).*randn(1,1) + a;   % initial drop size
%     
%     D(:,:,i) = droplet(height,drop_dim);     % simulate a water drop (size,???)
% end

%% Make empty vector for test points

test_H_mean = zeros(time,1); 
test_H = zeros(time,1);

%% Init. graphics
[surfplot,top] = initgraphics(Nvar);

%% Init. timer
tic;

%% Create ensamble of zeros/ones to store models

H = ones(Nvar+2,Nvar+2,Nens);   U = zeros(Nvar+2,Nvar+2,Nens);  V  = zeros(Nvar+2,Nvar+2,Nens);
Hx  = zeros(Nvar+1,Nvar+1,Nens); Ux  = zeros(Nvar+1,Nvar+1,Nens); Vx  = zeros(Nvar+1,Nvar+1,Nens);
Hy  = zeros(Nvar+1,Nvar+1,Nens); Uy  = zeros(Nvar+1,Nvar+1,Nens); Vy  = zeros(Nvar+1,Nvar+1,Nens);


%% Run Shallow Water Model
for itime = 1 : time
    
    % Output first 5 loops
    if mod(itime, 5) == 0
        fprintf('Current run: %d \n',itime)
        % Increase loop message frequency
    elseif itime <= 10
        fprintf('Current run: %d \n',itime)
    end
    
    % Nested loop for Ensembles
    for k = 1 : Nens
        
         % initialize water drop
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
        


        %% Update plot
        if mod(k,Nens) == 0
            
            test_H_mean(itime)  = mean(H(16,16,:));
            test_H(itime)  = H(16,16,Nens);

            C = abs(U(i,j,k)) + abs(V(i,j,k));  % Color shows momemtum
            set(surfplot,'zdata',H(i,j,k),'cdata',C);
            set(top,'string',sprintf('step = %d',itime))
            drawnow
        end
        
          %% Calc Error
        for q = 1:Nvar
            for p = 1:Nvar
                for s = 1 : Nens
                    temp = temp +  (H(q+1,p+1,s)-ObsValuesH(itime,q,p))^2;
                end
                RMS_H(itime,q,p) = sqrt(temp/Nens);
                temp = 0;
            end
        end
    end    
end
disp('Run time.....');
toc;

%% Save result for plotting and post-analysis purpose
filename = 'SWM_standard.mat';
save (filename);

xtime = 1:time;

for i = 1 : time
    RMSE(i) =  mean(mean(RMS_H(i,:,:)));
end


%% Plot results and error
figure(1)
plot(xtime,RMSE,'r',xtime,RMSE_EnKF,'b','LineWidth',1)
%axis([0 time -0.2 0.5])
title('mean of RMS Error of Reference Model')
xlabel('time')
ylabel('RMSE')

Obs_point = ObsValuesH(xtime,16,16);

figure(2)
plot(xtime,test_H_mean,'b',xtime,Obs_point,'r') %,'LineWidth',1
legend('Mean Height','Observed Height')
legend
% axis([0 time 0 5])
title('STANDARD Compare ObsValuesH to H_mean')
xlabel('time')
ylabel('value')

figure(3)
plot(xtime,test_H,'b',xtime,Obs_point,'r') % ,'LineWidth',1
legend('Rand Ens Height','Observed Height')
axis([0 time 0 5])
title('STANDARD Compare ObsValuesH to some ens of H')
xlabel('time')
ylabel('value')
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