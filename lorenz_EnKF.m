%
% cycle the ensemble Kalman filter on the Lorenz attractor.
%
% options for:
% - perturbed observations
% - "Potter" update (square too; all at once)
% - EnKF square-root alrgorithm (sequential obs with localization and inflation)
% - time-averaged observations
%
% Originator: G.J. Hakim, University of Washington
%             hakim@atmos.washington.edu
%
% Bug fixes: Ryan Torn, University at Albany, SUNY
%
% released under GNU General Public License version 3. http://www.gnu.org/licenses/gpl.html
%
% version control:
% $Date: 2012-11-21 14:14:08 -0800 (Wed, 21 Nov 2012) $
% $Revision: 125 $
% $Author: hakim $
% $Id: L63_ensda.m 125 2012-11-21 22:14:08Z hakim $

clear all
randn('state',0); % insure the same sequence of random numbers

update = 1;     % DA method (1= perturbed obs; 2= potter, 3=EnKF)
tavg = 0;       % time averaged observations (0 = no; 1 = yes)
q = .001;       % model error variance (covariance model is white for now)
R = eye(3)*1e-2;% observation error covariance
nassim = 10;   % number of assim cycles
ntimes = .1;    % do assimilation every ntimes nondimensional time units
climo = 0;      % option for climo covariance (0 = flow dependent)
Nens = 100;      % number of ensemble members
inflate = 1.1;  % covariance inflation 
H = eye(3);     % obs operator (eye(3) gives identity obs)
ind = 1;        % state vector index for archiving states
tau = 2;        % forecast lead time in units of ntimes (set 0 to turn off)

% control parameter settings for L63
par = [10 28 8/3];

% test parameter settings
testpar = [10 27 8/3];

if update == 1
  disp('running perturbed observations...')
  inflate = 0;
elseif update == 2
  disp('running the potter algorithm...')
elseif update == 3
  disp('running the enkf algorithm...')
else
  error('wrong update code specified')
end
  
% get a state on the attractor
disp('running onto the attractor...')
x = [10 20 30]; % initial conditions
[tims states] = ode45(@derivsL63,[0 100],x,[],par);

% IC for truth taken from last time (column vector):
xt = states(end,:)'; 

% populate initial ensemble analysis by perturbing true state
[tmp, Xa] = meshgrid(ones(1,Nens),xt); % 3 x Nens ensemble state matrix
pert = 1e-1*randn(3,Nens); 
Xa = Xa + pert; 

% check to see if we're using a fixed B
if climo > 0
  load L63_climo_B; Bc = B; 
end

disp(['Cycling ON the attractor...']);

% initialize history variables
hist_ver = []; hist_xbm = []; hist_xam = []; Xa_save = []; Xaf = [];
for k=1:nassim

  disp(['assimilation step ' int2str(k)])

  % step to the next assimilation time
  Xb = Xa; 
  
  % advance truth with the full nonlinear model
  [tims states] = ode45(@derivsL63,[0 ntimes],xt,[],par); xt = states(end,:)'; 

  % new observations from noise about truth; set verification values
  if tavg == 0
	 Y = H*xt + diag(diag(randn(3,1))*sqrtm(R));
	 ver = xt;
  else
	 Y = H*mean(states,1)' + diag(diag(randn(3,1))*sqrtm(R));
	 ver = mean(states,1)';
  end
  Nobs = size(Y,1);

  % advance background ensemble with the full nonlinear model
  for n = 1:1:Nens
	 [tims states] = ode45(@derivsL63,[0 ntimes],Xb(:,n),[],par); 
	 if tavg == 0
		Xb(:,n) = states(end,:)'; 
	 else
		Xb(:,n) = mean(states,1)'; % time mean ensemble
		Xbtpert(:,n) = states(end,:)' - Xb(:,n); % perts from time mean at final t
	 end		
  end
  
  xbm = mean(Xb,2);  % ensemble mean

  % new background error covariance from ensemble estimate
  if climo == 0
	 % remove the ensemble mean (XbM is a matrix version of xbm)
	 [tmp Xbm] = meshgrid(ones(Nens,1),xbm);
	 Xp = Xb - Xbm; 
	 if inflate ~= 0 % square-root filter
		Xbp = inflate*Xp;
		% additive zero-mean white model error
		Xbp = Xbp + q*randn(3,Nens);
	 else % perturbed obs
		Xbp = Xp;
	 end
	 B = (Xbp*Xbp')/(Nens - 1) + q*eye(3);
  else
	 B = Bc;
  end
  
  % update step for ensemble mean (xam), perturbations from the mean (Xap), 
  % and the full state Xa = Xam + Xap
  if update == 1
	 Xa = perturbed_obs(Xb,B,Y,R); 
	 xam = mean(Xa,2);
	 [tmp Xam] = meshgrid(ones(Nens,1),xam);
	 Xap = Xa - Xam; 
  elseif update == 2
	 [xam,Xap] = potter(xbm,Xbp,Y,H,R); 
	 [tmp Xam] = meshgrid(ones(Nens,1),xam);
	 Xa = Xam + Xap;
  elseif update == 3
	 loc = ones(Nobs,3); % this does no localization; size: Nobs x dof
	 [xam,Xap] = enkf_update(xbm,Xbp,Y,H,R,loc);	 
	 [tmp Xam] = meshgrid(ones(Nens,1),xam);
	 Xa = Xam + Xap;
  end

  A = (Xap*Xap')/(Nens - 1);

  disp(['trace of B and A: ' num2str(trace(B),2) '    ' num2str(trace(A),2)])

  %error statistics for the ensemble mean
  xbe(:,k) = xbm - ver; 
  xae(:,k) = xam - ver; 
  xye(:,k) = Y - ver;
  xaev(:,k) = diag(A);
  xbev(:,k) = diag(B);

  % check for filter divergence
  if abs(xae(1,k)) > 10 & abs(xae(2,k)) > 10
	 error('filter divergence!')
  end
  
  % history (for plotting)
  hist_ver(end+1) = ver(ind);
  hist_xbm(end+1) = xbm(ind);
  hist_xam(end+1) = xam(ind);
  
  % if time-averaged, add back the perts from the time mean
  if tavg > 0
	 Xa = Xa + Xbtpert;
  end

  % option to make and evaluate forecasts
  if tau >0 
	 if size(Xa_save,1) == tau
		Xaf = squeeze(Xa_save(tau,:,:)); % past analysis used to forecast current time is oldest saved
	 end
	 % update old analyses (first index has current; move older times to larger indices)
	 for m = 2:1:tau
		if m <= size(Xa_save,1)+1
		  Xa_save(m,:,:) = Xa_save(m-1,:,:);
		end
	 end
	 Xa_save(1,:,:) = Xa;
	 if isempty(Xaf) == 0 
		for n = 1:1:Nens
		  [tims states] = ode45(@derivsL63,[0 tau],Xaf(:,n),[],par); 
		  Xf(:,n) = states(end,:)'; 
		end
		% evaluate forecast error with current analysis
		xfm = mean(Xf,2);  % ensemble mean
		xfe(:,k) = xfm - ver; % error in ensemble mean
	 end
  end % tau
end % end of assimilation loop

% ensemble-mean error variance for the background and analysis
xbe_var = var(xbe,0,2);
xae_var = var(xae,0,2);

% plot error statistics
lab(1) = 'x'; lab(2) = 'y'; lab(3) = 'z';

% absolute error in the ensemble mean for each dimension
figure(1); clf
for k = 1:1:3
  subplot(3,1,k)
  plot(abs(xbe(k,:)),'b-'); hold on
  plot(abs(xae(k,:)),'r-')
  c = ['mean background: ' num2str(mean(abs(xbe(k,2:end))),3) ' +/- ' num2str(std(abs(xbe(k,2:end))),3) '    mean analysis:   ' num2str(mean(abs(xae(k,2:end))),3) ' +/- ' num2str(std(abs(xae(k,2:end))),3)];
  yl = get(gca,'ylim');
  yoff = yl(2)*1.1;
  text(0,yoff,c);
  ylabel(lab(k),'fontweight','bold','fontsize',12)
  if k == 1
	 h = title('absolute error in mean'); set(h,'position',[99,12,1],'fontweight','bold','fontsize',12)
  elseif k == 3
	 xlabel('Assimilation Step','fontweight','bold','fontsize',12)
  end
end

% background and analysis error variance for first variable (x)
figure(2); clf
plot(xbev(1,3:end),'b-'); hold on
plot(xaev(1,3:end),'r-')
c = ['mean background: ' num2str(mean(xbev(3:end)),3) ' +/- ' num2str(std(xbev(3:end)),3) '    mean analysis:   ' num2str(mean(xaev(3:end)),3) ' +/- ' num2str(std(xaev(3:end)),3)];
yl = get(gca,'ylim');
yoff = yl(1) - yl(2)*.12;
text(0,yoff,c);
xlabel('Assimilation Step','fontweight','bold','fontsize',12)
ylabel('Error Variance','fontweight','bold','fontsize',12)
title('Ensemble Kalman Filter Error Variance (x)','fontweight','bold','fontsize',12)
set(gca,'linewidth',2,'fontweight','bold','fontsize',10)
legend('background','analysis')

% time traces of the history variables
figure(3); clf
plot(hist_ver,'k-'); hold on
plot(hist_xbm,'b-')
plot(hist_xam,'r-')

%save L63_ensda xbe xae 
