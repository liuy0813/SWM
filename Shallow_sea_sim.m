function [states] = waterwave () %Drop_height, time

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

% Parameters

% define conditions of ensamble
num_runs = 10;

n = 64;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           % plot interval
%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval

drop_dim = 21;
D = zeros(21,21,num_runs);  % create empty array for different drops
for i = 1 : num_runs
    a = 1;                  % min size
    b = 8;                  % max size
    height = (b-a).*rand(1,1) + a;   % initial drop size
    D(:,:,i) = droplet(height,drop_dim);     % simulate a water drop (size,???)
end
% Initialize graphics

% [surfplot,top] = initgraphics(n);


% Outer loop, restarts.
max = 500; % total time
sample = max/10; % max/n where n is desired number of samples
nstep = 0;
states = zeros(10,4,num_runs);
test_num = 1;
loc = [16,16]; % grid is 64X64


% Initialize arrays to store states
% Uarray = zeros(max,1);
% Varray = zeros(max,1);
% Harray = zeros(max,1);



while nstep < max
    
    % Create ensamble of zeros here
    
    
    H = ones(n+2,n+2,num_runs);   U = zeros(n+2,n+2,num_runs);  V  = zeros(n+2,n+2,num_runs);
    Hx  = zeros(n+1,n+1,num_runs); Ux  = zeros(n+1,n+1,num_runs); Vx  = zeros(n+1,n+1,num_runs);
    Hy  = zeros(n+1,n+1,num_runs); Uy  = zeros(n+1,n+1,num_runs); Vy  = zeros(n+1,n+1,num_runs);
    
    
    %ndrop = ceil(rand*ndrops);
    nstep = 0;
    
    % Inner loop, time steps.
    
    while nstep < max
        nstep = nstep + 1;
        
        % Debugging code!!!
        if mod(nstep, 20) == 0
            fprintf('Current run: %d \n',nstep)
        elseif nstep <= 5
            fprintf('Current run: %d \n',nstep)
        end
        % initialize water drops
        for k = 1 : num_runs
            if nstep == 1;
                w = size(D(:,:,k),1);
                i = 5 +(1:w);
                j = 5 +(1:w);
                H(i,j,k) = H(i,j,k) + 0.5*D(:,:,k);
            end
            
            % Reflective boundary conditions
            H(:,1,k) = H(:,2,k);      U(:,1,k) = U(:,2,k);       V(:,1,k) = -V(:,2,k);
            H(:,n+2,k) = H(:,n+1,k);  U(:,n+2,k) = U(:,n+1,k);   V(:,n+2,k) = -V(:,n+1,k);
            H(1,:,k) = H(2,:,k);      U(1,:,k) = -U(2,:,k);      V(1,:,k) = V(2,:,k);
            H(n+2,:,k) = H(n+1,:,k);  U(n+2,:,k) = -U(n+1,:,k);  V(n+2,:,k) = V(n+1,:,k);
            
            %% Take a half time step to estimate derivatives at middle time.
            
            % x direction
            i = 1:n+1;
            j = 1:n;
            
            % height
            Hx(i,j,k) = (H(i+1,j+1,k)+H(i,j+1,k))/2 - dt/(2*dx)*(U(i+1,j+1,k)-U(i,j+1,k));
            
            % x momentum
            Ux(i,j,k) = (U(i+1,j+1,k)+U(i,j+1,k))/2 -  ...
                dt/(2*dx)*((U(i+1,j+1,k).^2./H(i+1,j+1,k) + g/2*H(i+1,j+1,k).^2) - ...
                (U(i,j+1,k).^2./H(i,j+1,k) + g/2*H(i,j+1,k).^2));
            
            % y momentum
            Vx(i,j,k) = (V(i+1,j+1,k)+V(i,j+1,k))/2 - ...
                dt/(2*dx)*((U(i+1,j+1,k).*V(i+1,j+1,k)./H(i+1,j+1,k)) - ...
                (U(i,j+1,k).*V(i,j+1,k)./H(i,j+1,k)));
            
            % y direction
            i = 1:n;
            j = 1:n+1;
            
            % height
            Hy(i,j,k) = (H(i+1,j+1,k)+H(i+1,j,k))/2 - dt/(2*dy)*(V(i+1,j+1,k)-V(i+1,j,k));
            
            % x momentum
            Uy(i,j,k) = (U(i+1,j+1,k)+U(i+1,j,k))/2 - ...
                dt/(2*dy)*((V(i+1,j+1,k).*U(i+1,j+1,k)./H(i+1,j+1,k)) - ...
                (V(i+1,j,k).*U(i+1,j,k)./H(i+1,j,k)));
            % y momentum
            Vy(i,j,k) = (V(i+1,j+1,k)+V(i+1,j,k))/2 - ...
                dt/(2*dy)*((V(i+1,j+1,k).^2./H(i+1,j+1,k) + g/2*H(i+1,j+1,k).^2) - ...
                (V(i+1,j,k).^2./H(i+1,j,k) + g/2*H(i+1,j,k).^2));
            
            %% Now take a full step that uses derivatives at middle point.
            
            i = 2:n+1;
            j = 2:n+1;
            
            % height
            H(i,j,k) = H(i,j,k) - (dt/dx)*(Ux(i,j-1,k)-Ux(i-1,j-1,k)) - ...
                (dt/dy)*(Vy(i-1,j,k)-Vy(i-1,j-1,k));
            % x momentum
            U(i,j,k) = U(i,j,k) - (dt/dx)*((Ux(i,j-1,k).^2./Hx(i,j-1,k) + g/2*Hx(i,j-1,k).^2) - ...
                (Ux(i-1,j-1,k).^2./Hx(i-1,j-1,k) + g/2*Hx(i-1,j-1,k).^2)) ...
                - (dt/dy)*((Vy(i-1,j,k).*Uy(i-1,j,k)./Hy(i-1,j,k)) - ...
                (Vy(i-1,j-1,k).*Uy(i-1,j-1,k)./Hy(i-1,j-1,k)));
            % y momentum
            V(i,j,k) = V(i,j,k) - (dt/dx)*((Ux(i,j-1,k).*Vx(i,j-1,k)./Hx(i,j-1,k)) - ...
                (Ux(i-1,j-1,k).*Vx(i-1,j-1,k)./Hx(i-1,j-1,k))) ...
                - (dt/dy)*((Vy(i-1,j,k).^2./Hy(i-1,j,k) + g/2*Hy(i-1,j,k).^2) - ...
                (Vy(i-1,j-1,k).^2./Hy(i-1,j-1,k) + g/2*Hy(i-1,j-1,k).^2));
            
            %Store H,U,V
            
            
            % This determines EnKF sample time
            if mod(nstep,sample) == 0
                
                states(nstep/sample,:,k) = [nstep,H(loc(1),loc(2),k),U(loc(1),loc(2),k),V(loc(1),loc(2),k)];
                
            end
        end
        
        if(mod(nstep,sample) == 0)
            %make avg matrix for state - all 3 variables for all x,y
            Tstart = tic
            ens_avg = get_ens_avg(H,U,V);
            obs = get_obs(nstep,0); %update to change error
            t1 = toc(Tstart)
            
            %generate errors
            ens_err = get_ens_err(H,U,V,ens_avg,nstep);
            obs_err = get_obs_err(0);
            t2 = toc(Tstart)
            
            %calc kalman gain
            k_gain = get_k_gain(ens_err,obs_err);
            t3 = toc(Tstart)
            %ofset ens_avg
            analysis_change = get_ana_chng(ens_avg,obs,k_gain);
            t4 = toc(Tstart)
            
            %update state
            update_ens(analysis_change,H,U,V);
            t2 = toc(Tstart)
            
        end
        
        
        %get observation matrix for current time - with error
        %Get analysis result
        %Either reinitialize all ensemble results from analysis result
        
        
        % Update plot
        %        if mod(nstep,nplotstep) == 0
        %           C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
        %           t = nstep*dt;
        %           tv = norm(C,'fro');
        %           set(surfplot,'zdata',H(i,j),'cdata',C);
        %           set(top,'string',sprintf('Drop Size = %6.2f',height)) % t = %6.2f, tv = %6.2f
        %           drawnow
        %        end
        %
        if all(all(isnan(H))), break, end  % Unstable, restart
    end
end

close(gcf)

return
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
colormap(hsv)
top = title('Shallow Sea Sim');

return
end

function ens_sum = get_ens_sum(H,U,V)
sum_mat = zeros(64,64,3);
ensemble = size(H,3);
for x = 1:64
    for y = 1:64
        for i=1:ensemble
            sum_mat(x,y,1) = sum_mat(x,y,1) + H(x,y,i);
            sum_mat(x,y,2) = sum_mat(x,y,2) + U(x,y,i);
            sum_mat(x,y,3) = sum_mat(x,y,3) + V(x,y,i);
        end
    end
end
ens_sum = sum_mat;
end

function avg = get_ens_avg(H,U,V)
ensemble = size(H,3);
avg = get_ens_sum(H,U,V) ./ ensemble;
end

function obs = get_obs(time,error)
ensemble_num = 100; %MAGIC NUMBER
sum = zeros(64,64,3);
file = importdata('REF_matrix.mat','-mat');

for x=1:64
    for y=1:64
        for data=1:3
            for run = 1:ensemble_num
                sum(x,y,data) = sum(x,y,data)+file(run,time,x,y,data);
            end
        end
    end
end

avg = sum ./ ensemble_num;
%do errror
obs = avg;
end

function error_ens = get_ens_err(H,U,V,ens_sum,T)
total_err = zeros(64,64,3);
ensemble = size(H,3);
for i=1:ensemble
    for x=1:64
        for y = 1:64
            total_err(x,y,1) = total_err(x,y,1) + (H(x,y,i) - ens_sum(x,y,1))...
                .*(H(x,y,1) - ens_sum(x,y,1))^T;
            
            total_err(x,y,2) = total_err(x,y,2) + (U(x,y,i) - ens_sum(x,y,2))...
                .*(H(x,y,i) - ens_sum(x,y,2))^T;
            
            total_err(x,y,3) = total_err(x,y,3) + (V(x,y,i) - ens_sum(x,y,3))...
                .*(H(x,y,i) - ens_sum(x,y,3))^T;
        end
    end
end
error_ens = total_err/(ensemble-1);
end

function obs_err = get_obs_err(input)
obs_err = eye(3); %we don't know what to do here
end

function k_gain = get_k_gain(ens_err, obs_err)
H = eye(3,1);
k_gain = zeros(64,64,3);
size(ens_err)
size(obs_err)
obs_err
size(H)
for i = 1:64
    for j = 1:64
        k_gain(i,j,:) = ens_err(i,j,:) .* (H.*ens_err(i,j,:)+obs_err)^-1; %what is H^t
    end
end
end

function analysis_change = get_ana_chng(ens_avg,obs,k_gain)

analysis_change =k_gain(obs - ens_avg);
end

function update_ens(ana_chng,H,U,V)
ensemble = size(H,3);
for x=1:64
    for y = 1:64;
        for i=1:ensemble
            H(x,y,i) = H(x,y,i) + ana_chng(x,y,1);
            U(x,y,i) = U(x,y,i) + ana_chng(x,y,2);
            V(x,y,i) = V(x,y,i) + ana_chng(x,y,3);
        end
    end
end
end

