function compare_error(ens_size, time, da_freq)
amp_array = [1.5,3]

if nargin == 0
    ens_size = 20
    time = 500
    da_freq = 50
    DA_amp = 1.5
end

n = size(amp_array,2);

for i = 1 : n
    
    DA_amp = amp_array(i);
    
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
    DA_dyn_sim(ens_size,time,da_freq, DA_amp);
    toc
    
    x = 1:time;
    stnd = importdata('Data/no_DA_error.mat');
    size(stnd)
    static = importdata('Data/DA_error.mat');
    dynam = importdata('Data/DA_dyn_error.mat');
    
    figure(1)
    plot(x,stnd,'g')
    hold on
    plot(x,static,'r')
    plot(x,dynam,'b')
    title('RMSE with amp factor of', 'fontsize', 20, 'fontweight', ...
        'bold');
    name=['Data/fig',num2str(i),'.fig'];
    saveas(gca,name);
    hold off
end
end