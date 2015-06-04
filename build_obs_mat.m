function build_obs_mats()

fprintf('Loading reference matrix...\n')

mat = importdata('REF_matrix.mat','-mat'); % default var name is ref_mat


[ens,time,Nvar,~,vals] = size(mat);  % read in matrix


fprintf('Bulding Obs matrix for H...\n')
Ens_H = mat(:,:,:,:,1);
mean_H = squeeze(mean(Ens_H,1));

fprintf('Bulding Obs matrix for U...\n')
Ens_U = mat(:,:,:,:,2);
mean_U = squeeze(mean(Ens_H,1));

fprintf('Bulding Obs matrix for V...\n')
Ens_V = mat(:,:,:,:,3);
mean_V = squeeze(mean(Ens_H,1));


save('OBS_matrix_H.mat','mean_H')
save('OBS_matrix_U.mat','mean_U')
save('OBS_matrix_V.mat','mean_V')
end




% save('REF_RMS_avg.mat','RMS_vector')

% Calculate RMS
% for i = 1 : time
%     for j = 1:Nvar
%         for k = 1:Nvar
%     REF_mean(i,j,k) = mean(mat(:,i,j,k,1));
%         end
%     end
% end
% sum = zeros(time,1);
% temp = 0;
% for i = 1 : time
%     for j = 1:Nvar
%         for k = 1:Nvar
%     temp = temp + RMS_vector(i,j,k);
%         end
%     end
%     sum(i) = temp/(Nvar^2);
%     temp = 0;
% end