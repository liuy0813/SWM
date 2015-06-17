function plot_mat(time, H_mat) %, U_mat, V_mat

%% Init. graphics
[surfplot,top] = initgraphics(64);

i = 1:64;
j = 1:64;

size( H_mat)
% size( U_mat)
% size( V_mat)
for k = 1 : time
%C = squeeze(abs(U_mat(k,i,j)) + abs(V_mat(k,i,j)));  % Color shows momemtum
% size(C)
% size(H_mat(k,i,j))
set(surfplot,'zdata',squeeze(H_mat(k,i,j))); %'cdata',C
set(top,'string',sprintf('step = %d',time))
drawnow
end

end

function [surfplot,top] = initgraphics(n)

% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

clf
shg
set(gcf,'numbertitle','off','name','Shallow_water')
x = (0:n-1)/(n-1);
surfplot = surf(x,x,ones(n,n),zeros(n,n));
grid off
axis([0 1 0 1 -1 3])
caxis([-1 1])
shading faceted
c = (1:64)'/64;
cyan = [c*0 c c];
colormap(cyan)
top = title('Shallow Sea Sim Ensemble');

return
end