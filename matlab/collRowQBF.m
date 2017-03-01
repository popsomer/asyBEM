% Construct (part of) a row of the collocation matrix.
% Input
%   i		- Row number
%   par		- The structure containing k,par, and so on
%   [j1		- Start column index, can not be negative for windows like (-0.1,0.2)]
%   [j2		- End column index]
%   [wind	- @(tau) The window function, default ones]
%   [collx	- Collocation point. If defined, we calculate part of the coupling matrix and assume no singularity]
%   [shift  - Shift with respect to colltau(i)]
% Output
%   row		- (part of) a row of the collocation matrix
function row = collRowQBF(i, par, j1, j2, wind, collx, shift)

if ~exist('j1','var'), j1 = 1; end
if ~exist('j2','var'), j2 = par.N; end

if ~exist('wind','var') || isempty(wind)
    if isfield(par, 'difase') % One can include phase information in difase which is interpolated and exploited here.
        wind = @(taut) exp(2i*pi*par.k*interp1(par.Tase,par.difase,taut, 'spline', 'extrap'));
    elseif isfield(par, 'phase') % One can include exact phase information in phase which is exploited here.
%         wind = @(taut) arrayfun(@(tautest) exp(1i*par.k*par.phase(tautest)), taut);
        wind = @(taut) exp(1i*par.k*par.phase(taut));
%         wind = @(taut) exp(-1i*par.k*phase(taut));
    else
        wind = @(taut) ones(size(taut));
    end
end

if par.dbf == 1
    if isfield(par, 'quadR')
        % Use a simple Riemann sum with the given number of points per collocation points interval to compute the integral. 
        % Using the middle sum avoids points of nonanalyticity.
        qbf_x = linspace(-1, 1, par.quadR*2+1);
        step = qbf_x(2)-qbf_x(1);
        qbf_x = qbf_x(1:end-1) + step/2;
%         qbf_w = [linspace(0, 1, par.quadR), linspace(1, 0, par.quadR)]; % /par.N already included
        qbf_w = (1-abs(qbf_x))/sum(1-abs(qbf_x));
        istep = round(1/step);
    else
        qbf_x = [  -1.00000000000000  -0.50000000000000                  0   0.50000000000000   1.00000000000000];
        qbf_w = [   0.01666666666667   0.26666666666667   0.43333333333333   0.26666666666667   0.01666666666667];
        % Step distance between quadrature pointsfor qbf
        step = 1/2; istep = 1/step;
    end
    A = -1;
    dbf = 1;
elseif par.dbf== 3 % Weights computed in iepack/basis/make_bfinfo_periodic_spline_equi.m using quad_sf with lo_d = sqrt(2)/16*[1 4 6 4 1]
    if isfield(par, 'quadR')
        error('Not implemented')
    end
    qbf_x = (-4:4)/2;
    qbf_w = [  0.000022045855379   0.009876543209876   0.084744268077601  0.238447971781305   0.333818342151675  ...
        0.238447971781305  0.084744268077601   0.009876543209877   0.000022045855379];
    step = 4/(length(qbf_w)-1);
    istep = round(1/step);
    
    A = -2;
    dbf = 3;
end
% The size of the row
Lj = j2-j1+1;

% if isfield(par, 'quadR')
% Total number of points in one row: there are istep additional points per element and 4+1 points for the first element,
% minus the istep already counted. Leave out D = 2 around the collocation point for the regular part.
Nj = Lj*istep + (length(qbf_w) - istep);
tau = (A+(j1-1):step:A+(j1-1)+(Nj-1)*step)/par.N;

if exist('collx','var') && ~isempty(collx)
    % if only SLP: kernelVals = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((repmat(collx,1,length(tau))-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
    kernelVals = kernel(par, tau, collx, wind);           
else
    tc = par.colltau(i);
    if exist('shift', 'var'),     tc = par.colltau(i) + shift; end
    % if only SLP: kernelVals = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(tau)))-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
    kernelVals = kernel(par, tau, tc, wind);
end
u = [par.t(end-dbf:end-1)-1, par.t, par.t(2:dbf+1)+1]; % the extended knots
row = zeros(1,Lj); % We could adjust A here, but that would not allow using parfor in the loop over i and might cause copying A.
% Perform Sweldens quadrature on the given kernelvalues; this gives NaN on singularities that are removed later.
for l=1:Lj
    if dbf == 1
        row(l) = qbf_w * (kernelVals((l-1)*istep+1:(l-1)*istep+length(qbf_w))).'/par.N;
    elseif (l == Lj) && (j2 == par.N) && (j1 == 1) % the last element in the row is wrapped because the first cubic spline is not symmetric
        row(1) = qbf_w * (kernelVals((l-1)*istep+1:(l-1)*istep+length(qbf_w))).'/par.N;
    elseif (l == Lj) % Using windows: this would need l=0 so recalculate for A+j1-2 and Lj=1
        tau = ( (A+j1-2):step:(A+j1-2+step*(length(qbf_w)-1) ) )/par.N;
        if exist('collx','var') && ~isempty(collx)
% if only SLP: kernelVals = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((repmat(collx,1,length(tau))-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
            kernelVals = kernel(par, tau, collx, wind);
        else
% if only SLP: kernelVals = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(tau)))-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
            kernelVals = kernel(par, tau, tc, wind);
        end
        row(1) = qbf_w * kernelVals.'/par.N;
    else % The problem is the first cubic spline is not centered around 0, so symmetry in the discretization matrix is lost.
        row(l+1) = qbf_w * (kernelVals((l-1)*istep+1:(l-1)*istep+length(qbf_w))).'/par.N;
    end
    if isfield(par, 'secondKind') && dbf == 1 && (~exist('collx','var') || isempty(collx)) % Add mass matrix
        k = l +j1 -1;
        row(l) = wind(par.colltau(k))*par.secondKind(1)*((tc -u(k))/(u(k+1) -u(k))*((tc >= u(k)) && (tc < u(k+1))) ...
            + (u(k+2) -tc)/(u(k+2) -u(k+1))*((tc >= u(k+1)) && (tc < u(k+2)))) + par.secondKind(2)*row(l);
        if k <= par.dbf
            row(l) = row(l) + wind(par.colltau(k))*par.secondKind(1)*(tc -u(k+par.N))/(u(k+par.N+1) -u(k+par.N))*((tc >= u(k+par.N)) && ...
                (tc < u(k+par.N+1))) + (u(k+par.N+2) -tc)/(u(k+par.N+2) -u(k+par.N+1))*((tc >= u(k+par.N+1)) && (tc < u(k+par.N+2)));
        end
    end
end
if exist('collx','var') && ~isempty(collx)
    return % A singularity is not possible when calculating the coupling matrix
end
% It would not be advisable to return when i is not in [j1,j2], as there could be windows close by i not in [j1,j2].
for j=-1:dbf % Removes NaN due to singularity
    idx = mod(i+j-1,par.N)+1;
    idxRow = mod(i+j-1,par.N)+2-j1;
    if (idxRow <= 0) || (idxRow > Lj)
        % Do nothing: this point does not fall in the [j1, j2] under interest. It could happen with periodic windows that j1=1 and
        % j2 < par.N and then negative j points are mapped to [par.N-2:par.N].
        continue
    end
    row(idxRow) = 0;
    tau1 = u(idx);
    tau2 = u(idx+dbf+1);
    if (tau1 <= tc) && (tau2 >= tc)
        d = 0;
    else
        d = min(abs(tc-tau1),abs(tc-tau2));
    end
    % Move the collocation point closer to the interval [tau1 tau2] if possible, using periodicity
    tct = tc; % == mod(tc,1)
    if d > 1/2
        if tc > tau2
            tct = tc-1; % Don't do mod(tc-1,1) here.
        else
            tct = tc+1;
        end
    end
    if (tau2-tau1) >= 1
        d = 0;
    elseif ( (tau1 <= tct) && (tau2 >= tct) && (tau1 < tau2) ) || ( ((tau1 <= tct) || (tau2 >= tct)) && (tau1 > tau2) )
        d = 0;
    else
        d = min(abs([tau1-tct, tau2-tct, tau1+1-tct,  tau2+1-tct,  tau1-1-tct,  tau2-1-tct]));
    end
    if d == 0
        % Calculate an element using hp-type quadrature on each (possibly singular but) smooth part of a piecewise smooth and singular integrand.
        sigma = 0.15;
        mu = 1;
        if isfield(par, 'quadR')
            n = 7 + par.quadR;
        else
            n = 7;
        end
        minsize = 1e-10;
        q = mod([par.corners tc]-tau1,1)+tau1;
        
        a = q((q >= tau1) & (q <= tau2));
        v = sort(u(idx:idx+1+dbf));
        for ai = 1:numel(a)
            if ~any(abs(v-a(ai)) < minsize/2^7)
                v = sort([v a(ai)]);
            end
        end
        for jj=1:length(v)-1
            a_j = v(jj);
            b_j = v(jj+1);
            D1 = min(abs([tc-a_j tc-a_j+1 tc-a_j-1]));
            D2 = min(abs([tc-b_j tc-b_j+1 tc-b_j-1]));
            if D1 <= D2
                [x,w] = quadHp(a_j, b_j, n, sigma, mu, minsize);
            else
                [x,w] = quadHp(b_j, a_j, n, sigma, mu, minsize);
                w = -w;
            end
            x = x';
            % Evaluate a periodic spline in the points x.
            z2 = zeros(size(x));
            if dbf == 1 % Using bspline_eval is correct but slower than inlining:
                for ii=1:length(x)
                    %   Evaluate the b-spline of degree 1 determined by the knots t at x.
                    z2(ii) = z2(ii) + (x(ii)-u(idx))/(u(idx+1)-u(idx))*((x(ii) >= u(idx)) && (x(ii) < u(idx+1))) ...
                        + (u(idx+2) -x(ii))/(u(idx+2)-u(idx+1))*((x(ii) >= u(idx+1)) && (x(ii) < u(idx+2)));
                    if idx <= 1 % Add the periodic parts for the basis functions at the end
                        z2(ii) = z2(ii) + (x(ii)-u(idx+par.N))/(u(idx+par.N+1)-u(idx+par.N))*((x(ii) >= u(idx+par.N)) && (x(ii) < u(idx+par.N+1))) ...
                            + (u(idx+par.N+2) -x(ii))/(u(idx+par.N+2)-u(idx+par.N+1))*((x(ii) >= u(idx+par.N+1)) && (x(ii) < u(idx+par.N+2)));
                    end
                end
            else
                for ii=1:length(x)
                    z2(ii) = z2(ii) + bspline_eval(u(idx:idx+par.dbf+1), par.dbf, x(ii));
                    % Add the periodic parts for the basis functions at the end.
                    if idx <= dbf
                        z2(ii) = z2(ii) + bspline_eval(u(idx+par.N:idx+par.N+par.dbf+1), par.dbf, x(ii));
                    end
                end
            end
            % if only SLP: f = wind(x).*1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(x)))-par.par(x) ).^2, 1)) ).*par.gradnorm(x).*z2;
            % If only SLP: row(idxRow) = row(idxRow) + f*w;
            row(idxRow) = row(idxRow) + (kernel(par, x, tc, wind).*z2)*w;
        end
    else % Use the general routine.
        if tau1 > tau2
            tau1 = tau1-1;
        end
        tau = (qbf_x+1)/2*(tau2-tau1)+tau1;
        % if only SLP: f = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(tau)))-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
        % if only SLP: row(idxRow) = qbf_w*f.'/par.N;
        row(idxRow) = qbf_w*kernel(par,tau,tc,wind).'/par.N;
    end
    
    if isfield(par, 'secondKind') && dbf == 1 && (~exist('collx','var') || isempty(collx)) % Add mass matrix
        k = idxRow +j1 -1;
        row(idxRow) = wind(par.colltau(k))*par.secondKind(1)*((tc -u(k))/(u(k+1) -u(k))*((tc >= u(k)) && (tc < u(k+1))) ...
            + (u(k+2) -tc)/(u(k+2) -u(k+1))*((tc >= u(k+1)) && (tc < u(k+2)))) + par.secondKind(2)*row(idxRow);
        if k <= par.dbf
            row(idxRow) = row(idxRow) + wind(par.colltau(k))*par.secondKind(1)*(tc -u(k+par.N))/(u(k+par.N+1) -u(k+par.N))*((tc >= u(k+par.N)) ...
                && (tc < u(k+par.N+1))) + (u(k+par.N+2) -tc)/(u(k+par.N+2) -u(k+par.N+1))*((tc >= u(k+par.N+1)) && (tc < u(k+par.N+2)));
        end
    end
    
end

end

% Evaluate the b-spline of degree k and determined by the knots t, at y.
function z = bspline_eval(t, k, y)

if k==0
    if (y >= t(1)) && (y < t(2))
        z = 1;
    else
        z = 0;
    end
else
    z = (y-t(1))/(t(end-1)-t(1)) * bspline_eval(t(1:end-1), k-1, y) ...
        + (t(end)-y)/(t(end)-t(2)) * bspline_eval(t(2:end), k-1, y);
end

end

% Evaluate the kernel for the single or double layer potential.
% This function was inlined for time measurements for integral formulations of the first kind.
function ker = kernel(par,tau,tc,wind)

if length(tc) == 1
    collx = par.par(tc*ones(size(tau)));
else
    collx = repmat(tc,1,length(tau));
end
if isfield(par, 'secondKind') % dlp
    xy = collx-par.par(tau);
    res = sqrt(sum(xy.^2, 1));
    z = 0*res;
    
    I1 = find(abs(res) > 1e-12); % I1 can be empty
    z(I1) = sum(par.normal(tau(I1)) .*xy(:,I1)) .*1i/4 .*par.k .*besselh(1, 1, par.k*res(I1)) ./res(I1);
    
    I2 = find(abs(res) <= 1e-12); % I2 can be empty
    g = par.grad(tau(I2));
    dg = par.derGrad(tau(I2));
    z(I2) = (dg(1,:).*g(2,:) -g(1,:).*dg(2,:))./(4*(g(1,:).^2 +g(2,:).^2).^(3/2) *pi);
    ker = wind(tau).*z.*par.gradnorm(tau);
else % slp
    ker = wind(tau).*1i/4.*besselh(0, 1, par.k*sqrt(sum((collx-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
end
end