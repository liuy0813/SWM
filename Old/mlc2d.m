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

dx=Lx/(nx-1);
dy=Ly/(ny-1);
% set up finite-difference mesh or grid:



[x y] = meshgrid(linspace(0,Lx,nx),linspace(0,Ly,ny));

whos

%plot(x,y, 'r.')  %plot shows resolution of Grid.


%%  Move thru time.  Make plots.
dt=.5;
Nsteps=100;
tests = 5;

% store states in arrays
leng = (ny-2)^2*Nsteps;
uarray = zeros(leng,tests);
varray = zeros(leng,tests);
harray = zeros(leng,tests);


for k = 1 : tests
    % initialize new run
    index = 1;
    %% Spaces for the etas.


eta=zeros(nx,ny);  u=zeros(nx,ny);   v=zeros(nx,ny);
etam=zeros(nx,ny);  um=zeros(nx,ny);  vm=zeros(nx,ny);  
etap=zeros(nx,ny); up=zeros(nx,ny);  vp=zeros(nx,ny);


%% reflective boundaries  
eta(:,1) = eta(:,2);      up(:,1) = up(:,2);       vp(:,1) = -vp(:,2);
eta(:,ny) = eta(:,ny-1);  up(:,ny) = up(:,ny-1);   vp(:,ny) = -vp(:,ny-1);
eta(1,:) = eta(2,:);      up(1,:) = -up(2,:);      vp(1,:) = vp(2,:);
eta(nx,:) = eta(nx-1,:);  up(nx,:) = -up(nx-1,:);  vp(nx,:) = vp(nx-1,:);




%% Initial condition. etam at t=0.

%Move the initial column of water around by changing io and jo. 
%k will change the width of the column
io=78;
jo=15;

% make random drop with magnitude 4<x<6
a = 4;
b = 6;
h = (b-a).*rand(1,1) + a; % Change height of initial drop here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


% simulates water drop to initialize
for i=1:nx 
   for j=1:ny
   
       etam(i,j)=10*exp((-((i-io)^2 + (j-jo)^2))/(h^2));
       
   end
end

eta=etam;
%this is what your initial condition looks like...
surf(x,y,eta)

for n=1:Nsteps
    t=n*dt;


    
    for i=2:nx-1
        for j=2:ny-1
    up(i,j)=um(i,j)-g*(dt/dx)*(eta(i+1,j)-eta(i-1,j));
    vp(i,j)=vm(i,j)-g*(dt/dy)*(eta(i,j+1)-eta(i,j-1));
    
    etap(i,j)=2*eta(i,j)-etam(i,j)+(((2*dt^2/dx*dy))*...
            (eta(i+1,j)+eta(i-1,j)+eta(i,j+1)+eta(i,j-1)-4*eta(i,j)));
        
        % Store state variables
        uarray(index,k) = up(i,j);
        varray(index,k) = vp(i,j);
        harray(index,k) = etap(i,j);
        
        index = index + 1;

        end
    end
    


 etam=eta;
 eta=etap;

pause(0.1)
figure(1)



subplot(1,2,1)
z=etap;
surf(x,y,eta)

colormap('Gray')
zlim([-5 5])
shading interp

zlabel('H +/- \eta')

subplot(1,2,2)
surf(x,y,eta)

view(0,90)
shading flat

xlabel('Width (m)')
ylabel('Width (m)')



end

end

% take sample states at every 1000th interval 
indexTwo = 1;
states = zeros((ceil(index/100)),3);
length(states)

statesU = fopen('statesU.txt', 'wt');
statesV = fopen('statesV.txt', 'wt');
statesH = fopen('statesH.txt', 'wt');


for i=1:(index/100)
 states(i,:) = [uarray(indexTwo),varray(indexTwo),harray(indexTwo)];
fprintf(statesU,'%f',uarray(indexTwo,:));
fprintf(statesU,'\n');
fprintf(statesV,'%f',varray(indexTwo,:));
fprintf(statesV,'\n');
fprintf(statesH,'%f',harray(indexTwo,:));
fprintf(statesH,'\n');
%  fprintf(fileID,formatSpec,A1,...,An)
%  fprintf(fileID,formatSpec,A1,...,An)
 indexTwo = indexTwo + 1000;
end






