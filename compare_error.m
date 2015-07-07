function compare_error(ens_size, time, da_freq)
if nargin == 0
    ens_size = 20
    time = 500
    da_freq = 20
end

tic
fprintf('Running static DA')
tic
DA_sim(ens_size,time,da_freq);
toc

fprintf('Running standard')
tic
SWM_sim(ens_size,time);
toc

fprintf('Running dynamic DA')
tic
DA_dyn_sim(ens_size,time,da_freq);
toc

x = 1:time;
stnd = importdata('Data/no_DA_error.mat');
size(stnd)
static = importdata('Data/DA_error.mat');
dynam = importdata('Data/DA_dyn_error.mat');


plot(x,stnd)
hold on
plot(x,static)
plot(x,dynam)
end