%% Fast, efficient approximation for SWE.


clear all
clc

%% This is my first cell, create constants, and grid.
%Constants
g=9.81; % units of gravity are m/s^2


%Grid length, number and spacing
Lx=100;
Ly=100;
nx=101;
ny=101;

dx=Lx/(nx-1); % dx & dy = 1
dy=Ly/(ny-1);

% set up finite-difference mesh or grid:

dt=.5;
Nsteps=10;
runs = 3;
uarr = zeros(ny,ny,Nsteps/dt,runs);
varr = zeros(ny,ny,Nsteps/dt,runs);
etaarr = zeros(ny,ny,Nsteps/dt,runs);
for k=1:runs

[x y] = meshgrid(linspace(0,Lx,nx),linspace(0,Ly,ny));

whos;

%plot(x,y, 'r.')  %plot shows resolution of Grid.

%% Spaces for the etas, i.e. sea-level height above or below mean sea-level location

h=zeros(nx,ny);  u=zeros(nx,ny);   v=zeros(nx,ny);
hx=zeros(nx,ny);  ux=zeros(nx,ny);  vx=zeros(nx,ny);
hy=zeros(nx,ny); uy=zeros(nx,ny);  vy=zeros(nx,ny);


%% reflective boundaries
h(:,1) = h(:,2);      uy(:,1) = uy(:,2);       vy(:,1) = -vy(:,2);
h(:,ny) = h(:,ny-1);  uy(:,ny) = uy(:,ny-1);   vy(:,ny) = -vy(:,ny-1);
h(1,:) = h(2,:);      uy(1,:) = -uy(2,:);      vy(1,:) = vy(2,:);
h(nx,:) = h(nx-1,:);  uy(nx,:) = -uy(nx-1,:);  vy(nx,:) = vy(nx-1,:);




%% Initial condition. etam at t=0.

%Move the initial column of water around by changing io and jo.
%k will change the width of the column
io=78;
jo=15;
h_start=rand()+5


% simulates water drop to initialize
for i=1:nx
   for j=1:ny

       hx(i,j)=10*exp((-((i-io)^2 + (j-jo)^2))/(h_start^2));

   end
end


h=hx;

%this is what your initial condition looks like...
surf(x,y,h)
%%  Move thru time.  Make plots.

% store states in arrays


%index
index = 1;

for n=1:Nsteps
    t=n*dt;



    for i=2:nx-1
        for j=2:ny-1
    uy(i,j)=ux(i,j)-g*(dt/dx)*(h(i+1,j)-h(i-1,j));
    vy(i,j)=vx(i,j)-g*(dt/dy)*(h(i,j+1)-h(i,j-1));

    hy(i,j)=2*h(i,j)-hx(i,j)+(((2*dt^2/dx*dy))*...
(h(i+1,j)+h(i-1,j)+h(i,j+1)+h(i,j-1)-4*h(i,j)));

        % Store state variables


        uarr(i,j,n,k) = uy(i,j);
        varr(i,j,n,k) = vy(i,j);
        etaarr(i,j,n,k) = hy(i,j);

        end
    end



 hx=h;
 h=hy;

pause(0.1)
figure(1)

subplot(1,2,1)
z=hy;
surf(x,y,h)

colormap('Gray')
zlim([-5 5])
shading interp

zlabel('H ± \eta')

subplot(1,2,2)
surf(x,y,h)

view(0,90)
shading flat

xlabel('Width (m)')
ylabel('Width (m)')



end

% take sample states at every 100th interval
indexTwo = 1;
states = zeros((ceil(index/100)),3);
length(states);



end
for l = 1 : runs
    etaarr(75,15,7,l)
end

