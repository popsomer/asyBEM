% Possibly sort, periodize and join the bounds and then construct a row of the collocation matrix using windows.
% Parts where the window == 1 can only become larger by joining windows. multElWind and multKernelWind are not mutually exclusive:
% if they are both zero, there is a block window.
% Input
%   i              - Row number, or its negative if the bounds do not have to be sorted nor joined any more
%   par            - The structure containing k,par and so on
%   bounds         - Each column gives the lambda, rho, l-lambda and rho-r of the window, but these could overlap or be too small.
%   [multElWind    - nonzero when multiplying the original matrix elements with the window function evaluated on the respective support, default 1]
%   [multKernelWind- nonzero when multiplying the kernel with the window function, default 0]
%   [collx         - Collocation point for when computing row in coupling matrix]
%   [normRois      - rois(i,:) divided by the threshold, in case windows are computed using correlation 'distances']
%   [ftz           - the factor below the threshold below which the window is zero, ftz=1 means block window, when correlation 'distances']
%   [roit          - the scaling between rois and A2, when correlation 'distances']
% Output
%   row            - a row of the collocation matrix
function row = windRow(i,par,bounds,multElWind, multKernelWind, collx,normRois, ftz,roit)
if ~exist('multElWind','var') || isempty(multElWind)
    multElWind = 1;
end
if ~exist('multKernelWind','var') || isempty(multKernelWind)
    multKernelWind = 0;
end
if ~exist('collx','var')
    collx = [];
end
windGreen = @(t) ones(size(t)); % The modification to the Green's function

%% Compute windows using correlation 'distances'
if exist('normRois','var') && ~isempty(normRois)
    getjr = @(t) arrayfun(@(tti) find(min(abs(roit-tti)) == abs(roit-tti),1,'first'), t);
    if multKernelWind
        windGreen = @(t) chi(abs(normRois(getjr(t))), ftz, 1,Inf,Inf,0);
    end
    if ~isempty(collx)
        row = collRowQBF('error', par, 1, par.N, windGreen, collx);
    else
        row = collRowQBF(i, par, 1, par.N, windGreen);
    end
    if multElWind
        row = chi(abs(normRois(getjr(par.colltau))), ftz, 1,Inf,Inf,0).*row;
    end
    return
end

%% Possibly sort the bounds
if i < 0
    i = -i;
else
    bounds = sortBounds(i,par,bounds, collx);
end

%% Construct the row
% Also take colltau into account: when [-1e-12; 0.9999;...] then all collocation points fall into the window, so we make the full row.
if isempty(bounds) || (max(bounds(2,:)-bounds(1,:)) > par.colltau(end)) || (norm(bounds) == 0)
    if exist('collx','var') && ~isempty(collx)
        row = collRowQBF('error',par,1,par.N,[],collx);
    else
        row = collRowQBF(i,par);
    end
    return
end
row = zeros(1,par.N);
for nrWind = 1:size(bounds,2)
    lambda = bounds(1,nrWind);
    rho = bounds(2,nrWind);
    l = lambda + bounds(3,nrWind);
    r = rho - bounds(4,nrWind);
    if lambda > rho
        error('Expected rho >= lambda');
    elseif l > r
        error('Constant 1 part of window wrong in bounds');
    elseif (nrWind >= 2) && (lambda < bounds(2,nrWind-1))
        error('Bounds were not correctly joined');
    end
    
    scp = @(t) chi(t,lambda,l, r,rho,1); % Smooth Cutoff Periodic
    % We could test whether scp(par.colltau) has min 0 and max 1 but this might not be the case for a small window or l close to r.
    if multKernelWind
        windGreen = scp;
    end
    [j1, j2, noSplit] = bounds2Ind(par.colltau,lambda,rho);
    if noSplit % j1 = j2 or j2-1 should have have been removed when windows are too small, so now colltau(j1) < lambda < rho < colltau(j2)
        if  exist('collx','var') && ~isempty(collx)
            row(j1:j2) = collRowQBF('error', par, j1, j2, windGreen, collx);
        else
            row(j1:j2) = collRowQBF(i, par, j1, j2, windGreen);
        end
        if multElWind
            row(j1:j2) = scp(par.colltau(j1:j2)).*row(j1:j2);
        end
    else % A window like (-0.1,0.2) has to be split.
        if  exist('collx','var') && ~isempty(collx)
            row(j1:par.N) = collRowQBF('error',par,j1,par.N,windGreen,collx);
            row(1:j2) = collRowQBF('error',par,1,j2,windGreen,collx);
        else
            row(j1:par.N) = collRowQBF(i,par,j1,par.N,windGreen);
            row(1:j2) = collRowQBF(i,par,1,j2,windGreen);
        end
        if multElWind
            row(j1:par.N) = scp(par.colltau(j1:par.N)).*row(j1:par.N);
            row(1:j2) = scp(par.colltau(1:j2)).*row(1:j2);
        end
    end
end

end
