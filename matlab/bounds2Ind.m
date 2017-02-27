% Convert the bounds to collocation indices.
% Input
%   ct          - Collocation points (par.colltau)
%   lambda, rho - Window goes from lambda to rho
% Output
%   j1,j2       - Required indices in ct are j1:j2 if noSplit, else 1:j2 and j1:end
%   noSplit     - 0 if there are two intervals due to periodicity
function [j1, j2, noSplit] = bounds2Ind(ct, lambda, rho)
j1 = find(ct <= lambda, 1, 'last');
if ((lambda < 0) && (lambda+1 > ct(end))) || ((lambda > ct(end)) && (lambda < 1))
	j1 = 1; % j1 = length(ct); gives noSplit = 0 so that would try j1:end while we need 1:j2
elseif lambda < 0
	j1 = find(ct <= lambda+1, 1, 'last');
elseif lambda > ct(end)
	j1 = find(ct <= lambda-1, 1, 'last'); % lambda-1 cannot fall between 0 and ct(1)=0
end
j2 = find(ct >= rho, 1, 'first');
if (rho < 0) && (rho+1 > ct(end))
	j2 = length(ct);
elseif rho < 0
	j2 = find(ct >= rho+1, 1, 'first');
elseif (rho > ct(end)) && rho < 1
	j2 = length(ct); % This was added because else j2=1 and noSplit =0 so that would try 1:1 instead of j1:N.
elseif rho > ct(end)
	j2 = find(ct >= rho-1, 1, 'first');
end
if isempty(j1) || isempty(j2)
	error(['bounds2Ind failed: j1 = ' num2str(j1) ', j2 = ' num2str(j2) ', lambda = ' num2str(lambda) ', rho = ' num2str(rho)]);
end
noSplit = (j1<j2) && ~(( (ct(j1) < rho) && (rho < 1+lambda) && (1+lambda < ct(j2)) ) || ...
	((ct(j1)+1 < rho) && (rho < 1+lambda) && (1+lambda < ct(j2)+1)));
