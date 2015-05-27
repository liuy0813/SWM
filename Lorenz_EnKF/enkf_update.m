
function [xam,Xap] = enkf_update(xbm, Xbp, Y, H, R, loc)
%Inputs
%  Xbp = ensemble of states
% xam = ensemble mean
% H identity function
% R model Error matrix
% Y is somthing! (prior guess)


%
% Originator: G.J. Hakim, University of Washington
%             hakim@atmos.washington.edu
%
% Modified by Ryan Torn, University at Albany, SUNY
%  
% released under GNU General Public License version 3. http://www.gnu.org/licenses/gpl.html
%
% version control:
% $Date: 2012-11-20 11:48:50 -0800 (Tue, 20 Nov 2012) $
% $Revision: 124 $
% $Author: hakim $
% $Id: enkf_update.m 124 2012-11-20 19:48:50Z hakim $

Nens = size(Xbp,2); % number of ensemble members
Ndim = size(Xbp,1); % number of degrees of freedom in state vector
Nobs = size(Y,1);  % number of observations

xam = xbm;
Xap = Xbp;

% loop over observations
for ob = 1:Nobs

  % vector of model estimated obs with mean removed
  Ye  = H(ob,:) * (repmat(xam, [1,Nens]) + Xap);
  mye = mean(Ye,2); 

  ye = Ye - mye;
  varye = var(ye);
  obs_err = R(ob,ob);

  % innovation
  innov = Y(ob) - mye;

  % innovation variance
  kdenom = (varye + obs_err);

  % B H^T ~ X (H X)^T ~ cov(Xm, ye)
  kcov = Xap*ye'/(Nens-1);

  % localize the gain
  kcov = kcov .* loc(ob,:)';
  
  % kalman gain
  kmat = (kcov ./ kdenom);

  % update the mean
  xam = xam + (kmat * innov);

  % update the ensemble members
  beta =  (1 ./ (1 + sqrt(obs_err ./ (varye + obs_err))));
  kmat = beta * kmat;
  Xap = Xap - (kmat * ye);

end
