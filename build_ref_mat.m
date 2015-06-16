function build_ref_mat(Nens,time)

ref_mat = zeros(Nens,time,64,64);

% timing

n = 64;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;

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
%% Init. graphics
[surfplot,top] = initgraphics(n);



% Outer loop, restarts.
max = time; % total time
sample = max/max; % max/n where n is desired number of samples
loc = [2:65,1:65]; % grid is 64X64


% Initialize arrays to store states
% Uarray = zeros(max,1);
% Varray = zeros(max,1);
% Harray = zeros(max,1);

%Start timer
Tstart = tic;
fprintf('building matrix... \n');


% Create ensamble of zeros here
H = ones(n+2,n+2,Nens);   U = zeros(n+2,n+2,Nens);  V  = zeros(n+2,n+2,Nens);
Hx  = zeros(n+1,n+1,Nens); Ux  = zeros(n+1,n+1,Nens); Vx  = zeros(n+1,n+1,Nens);
Hy  = zeros(n+1,n+1,Nens); Uy  = zeros(n+1,n+1,Nens); Vy  = zeros(n+1,n+1,Nens);

% Matrix for storing vals
mean_H = zeros(time,64,64);

% Inner loop, time steps.

for nstep = 1 : time
    
    
    % Output progress message every 100 steps
    if mod(nstep, 100) == 0
        fprintf('Number of step: %d   Time: %2.5f \n', nstep, toc(Tstart))
    elseif (nstep <= 5)
        fprintf('Number of step: %d   Time: %2.5f \n', nstep, toc(Tstart))
    end
    
    % initialize water drops
    for k = 1 : Nens
        if nstep == 1;
            w = size(D(:,:,k),1);
            i = 5 +(1:w);
            j = 5 +(1:w);
            H(i,j,k) = H(i,j,k) + D(:,:,k);
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
        
        
        % Save ref_mat
%         if mode(nstep,5) == 0
        ref_mat(k,nstep,:,:) = H(i,j,k);
%         ref_mat(k,nstep,:,:,2) = U(i,j,k);
%         ref_mat(k,nstep,:,:,3) = V(i,j,k);
%         end
        
        % take average of ensemble
       
        
    end
     % take average of ensemble
    Ens_H = ref_mat(:,nstep,:,:);
    mean_H(nstep,:,:) = squeeze(mean(Ens_H,1));
    
    % Update plot
        C = abs(U(i,j,k)) + abs(V(i,j,k));  % Color shows momemtum
        set(surfplot,'zdata',H(i,j,k),'cdata',C);
        set(top,'string',sprintf('step = %d',nstep))
        drawnow

    
end

% % stop timer
% runtime = toc(Tstart);
% fprintf('Total time: %2.5f', runtime);

% save matrix
save('OBS_matrix_H.mat','mean_H')
% save('REF_matrix.mat','ref_mat')
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

% ------------------------------------

