function [states] = waterwave ( )

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
%  Reference:
%
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www

% Parameters

n = 64;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           % plot interval
%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval
a = 1.4;                  % min size
b = 1.6;                  % max size
height = (b-a).*rand(1,1) + a;   % initial drop size
D = droplet(height,21);     % simulate a water drop (size,???)

% Initialize graphics

[surfplot,top] = initgraphics(n);


% Outer loop, restarts.
max = 1000; % total time
sample = max/max; % max/n where n is desired number of samples
nstep = 0;
states = zeros(max/sample,4);
test_num = 1;
loc = [16,16]; % grid is 64X64


% Initialize arrays to store states
% Uarray = zeros(max,1);
% Varray = zeros(max,1);
% Harray = zeros(max,1);

while nstep < max
   
   H = ones(n+2,n+2);   U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
   Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
   Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
   %ndrop = ceil(rand*ndrops);
     nstep = 0;

   % Inner loop, time steps.

   while nstep < max
       nstep = nstep + 1;

       % initialize water drops
       if nstep < 2;
           w = size(D,1);
           i = 5 +(1:w);
           j = 5 +(1:w);
           H(i,j) = H(i,j) + 0.5*D;
       end
     
       % Reflective boundary conditions
       H(:,1) = H(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
       H(:,n+2) = H(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
       H(1,:) = H(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
       H(n+2,:) = H(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:);

       %% Take a half time step to estimate derivatives at middle time.
   
       % x direction
       i = 1:n+1;
       j = 1:n;
   
       % height
       Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
   
       % x momentum
       Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
                 dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                            (U(i,j+1).^2./H(i,j+1) + g/2*H(i,j+1).^2));
   
       % y momentum
       Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
                 dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
                            (U(i,j+1).*V(i,j+1)./H(i,j+1)));     
                        
       % y direction
       i = 1:n;
       j = 1:n+1;
   
       % height
       Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
   
       % x momentum
       Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
                            (V(i+1,j).*U(i+1,j)./H(i+1,j)));
       % y momentum
       Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
                 dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                            (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2));
   
       %% Now take a full step that uses derivatives at middle point.

       i = 2:n+1;
       j = 2:n+1;
   
       % height
       H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                         (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
       % x momentum
       U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
                         (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
                       - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
                         (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
       % y momentum
       V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
                         (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
                       - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
                         (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2));
   
       %Store H,U,V  
       if mod(nstep,sample) == 0
           states(test_num,:) = [nstep,H(loc(1),loc(2)),U(loc(1),loc(2)),V(loc(1),loc(2))];
           test_num = test_num+1;
           
       end

       % Update plot
       if mod(nstep,nplotstep) == 0
          C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
          t = nstep*dt;
          tv = norm(C,'fro');   
          set(surfplot,'zdata',H(i,j),'cdata',C);
          set(top,'string',sprintf('Drop Size = %6.2f',height)) % t = %6.2f, tv = %6.2f
          drawnow
       end
      
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
   colormap(winter)
   top = title('Shallow Sea Sim');

  return
end
