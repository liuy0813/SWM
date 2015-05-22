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



function P_f = error_forecast(ens_state, obs_state, time)
%change 1 to 2 if array is formatted differently
l = size(ens_state,1);
sum = 0;
for i = 1: l
    %do we really raise to the power of time????
    sum = sum + (ens_state(i)-obs_state)*(ens_state(i)-obs_state)^time;
end
avg = sum / l;
P_f = avg;
end

function k_gain = get_kalman(f_err, model_err)
%should model_err be hard coded?
k_gain = f_err.*(f_err + model_err)^-1;
end



