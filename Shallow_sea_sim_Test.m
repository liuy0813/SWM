% Varience test for Shallow Sea Simulation
% check specific spot on maps for STATE = (h[t,x,y],u[t,x,y],v[t,x,y])
% function Shallow_sea_sim()  set up take 10 samples
% n = number of tests
% states = 3D matrix of state results 
% x = test number
% y = h,u,v
% z = # simulation run

n = 100;
states = zeros(10000,4,n);

 for i=1:n
     states(:,:,i) =  Shallow_sea_sim();
 end
 
save('OBS_matrix.mat','states')