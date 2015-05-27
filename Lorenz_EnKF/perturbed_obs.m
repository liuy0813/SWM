function [Xa] = perturbed_obs(Xb,B,Y,R)

%
% Originator: G.J. Hakim, University of Washington
%             hakim@atmos.washington.edu
%
% released under GNU General Public License version 3. http://www.gnu.org/licenses/gpl.html
%
% version control:
% $Date: 2009-02-25 16:04:40 -0800 (Wed, 25 Feb 2009) $
% $Revision: 44 $
% $Author: hakim $
% $Id: perturbed_obs.m 44 2009-02-26 00:04:40Z hakim $

% update an ensemble by perturbed observations
  Nens = size(Xb,2); % number of ensemble members
  Ndim = size(Xb,1); % number of degrees of freedom in state vector

  for n = 1:1:Nens
	 xb = Xb(:,n);
	 Yn = Y + diag(diag(randn(Ndim,1))*sqrtm(R));
	 K = B*inv(B + R);
	 Xa(:,n) = Xb(:,n) + K*(Yn - Xb(:,n));
  end

