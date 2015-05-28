function ref_mat = build_ref_mat(ens,time)

ref_mat = zeros(ens,time,64,64,3);


n = 64;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.02;               % hardwired timestep
dx = 1.0;
dy = 1.0; 
%ndrops = 5;              % maximum number of drops
%dropstep = 500;          % drop interval

drop_dim = 21;
D = zeros(21,21,ens);  % create empty array for different drops
for i = 1 : ens
a = 1;                  % min size
b = 8;                  % max size
height = (b-a).*rand(1,1) + a;   % initial drop size
D(:,:,i) = droplet(height,drop_dim);     % simulate a water drop (size,???)
end
% Initialize graphics



% Outer loop, restarts.
max = time; % total time
sample = max/max; % max/n where n is desired number of samples
nstep = 0;
loc = [1:64,1:64]; % grid is 64X64


% Initialize arrays to store states
% Uarray = zeros(max,1);
% Varray = zeros(max,1);
% Harray = zeros(max,1);



while nstep < max
    
   % Create ensamble of zeros here
   

   H = ones(n+2,n+2,ens);   U = zeros(n+2,n+2,ens);  V  = zeros(n+2,n+2,ens);
   Hx  = zeros(n+1,n+1,ens); Ux  = zeros(n+1,n+1,ens); Vx  = zeros(n+1,n+1,ens);
   Hy  = zeros(n+1,n+1,ens); Uy  = zeros(n+1,n+1,ens); Vy  = zeros(n+1,n+1,ens);
  
   
   %ndrop = ceil(rand*ndrops);
     nstep = 0;

   % Inner loop, time steps.

   while nstep < max
       nstep = nstep + 1;
       
       % Debugging code!!!
       if mod(nstep, 100) == 0
           nstep
       end
       % initialize water drops
       for k = 1 : ens
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

           ref_mat(k,nstep,:,:,1) = H(loc(1),loc(2),k);
             ref_mat(k,nstep,:,:,2) = U(loc(1),loc(2),k);
             ref_mat(k,nstep,:,:,3) = V(loc(1),loc(2),k);
           
           
       end
       
        

%       
       if all(all(isnan(H))), break, end  % Unstable, restart
       end
  end
end
save('REF_matrix.mat','ref_mat')
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

% n = 100;
% states = zeros(10000,4,n);
% 
%  for i=1:n
%      states(:,:,i) =  Shallow_sea_sim();
%  end
%  
% save('OBS_matrix.mat','states')