function [ states ] = SWM_ekf( input_args )

% make empty array of length 20 for enssamble runs
drops = zeros(20,1);
time = [1000,2000,3000,4000,5000];

n = length(drops);
m = length(time);

for i = 1 : 20
drops(i) = (1.6-1.4).*rand(1,1) + 1.4;
end

states = zeros(10,4,n,m);

 for i=1:m
     for j = n
     states(:,:,j,i) =  Shallow_sea_sim(drops(j),time(i));
     end
 end



end

