function [ d ] = KLD( h1,h2 )
%# you may want to do some input testing, such as whether h1 and h2 are
%# of the same size

%# preassign the output
d = zeros(size(h1));

%# create an index of the "good" data points
goodIdx = h1>0 & h2>0; %# bin counts <0 are not good, either

d1 = sum(h1(goodIdx) .* log(h1(goodIdx) ./ h2(goodIdx)));
d2 = sum(h2(goodIdx) .* log(h2(goodIdx) ./ h1(goodIdx)));

%# overwrite d only where we have actual data
%# the rest remains zero
d(goodIdx) = d1 + d2;
end