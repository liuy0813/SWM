function [xam,Xap] = potter(xbm,Xbp,Y,H,R)

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
% $Id: potter.m 124 2012-11-20 19:48:50Z hakim $

% update an ensemble by the Potter method
  Nens = size(Xbp,2); % number of ensemble members
  Ndim = size(Xbp,1); % number of degrees of freedom in state vector
  Nobs = length(Y);

  xam = xbm;
  Xap = Xbp;

  % square-root update step by Potter algorithm
  B = (Xbp*Xbp')/(Nens - 1);
  for ob = 1:Nobs
	 F = (H(ob,:)*Xap)';
	 alph = 1.0/(F'*F./(Nens-1) + R(ob,ob));
	 gam = 1.0/(1.0+sqrt(alph*R(ob,ob)));
	 K = alph*Xap*F ./ (Nens-1);
	 xam = xam + K*(Y(ob) - (H(ob,:)*xam)); % update state
	 Xap = Xap - gam*K*F'; % update square root
%         disp([int2str(ob) ' ' num2str(Y(ob)) ' ' num2str(xbm(1) + Xbp(1)) ' ' num2str(xam(1) + Xap(1))]);
  end
  

