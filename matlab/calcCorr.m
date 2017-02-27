% Calculate the correlations, or the resulting bounds.
% Input
%   par       - The structure containing k,par, and so on
%   c         - The solution vector to compute the correlations with
%   Tcor      - The width of the windows
%   pd        - The percentage of support of the windows with a decay
%   [tim      - If present and nonempty, tim(1) is the number of seconds after which 
%                  to print tim(2) + expected days for this method, assuming tic has been started just before calling calcCorr]
%   [A        - If present, the system matrix to compute the correlations with]
%   [method   - 'colR' if R computed columnwise from A (default), 'colB' if bounds computed columnwise from A (faster than rowB), 
%               'rowR' if R computed rowwise from A, 'rowB' if bounds computed rowwise from A, 'colRint' if A unused,...]
%   [a2w      - For \tilde{A} around the Green singularity: w(t_i,tau) = chi(tau-t_i, -a2w(1)/k, (a2w(2)-1)*a2w(1)/k, (1-a2w(2))*a2w(1)/k, 
%               a2w(1)/k) and a2w(2)*a2w(1)/k is the width of the decaying part left and right of where corr \approx threshold.
% Output
%   R         - The correlations where each row corresponds to par.colltau
%   sigma     - The locations of the windows, corresponding the columns of R
%   [obbounds - Additional information for multiple scattering obstacles]
function [R, sigma, obbounds] = calcCorr(par, c, Tcor, pd, tim, A, method, a2w)
if ~exist('tim','var') || isempty(tim)
    tim = [inf, 0];
end
if ~exist('method','var')
    method = 'colR';
end

fr = 1.5; % factor of number of columns over rows of R
% fr = 1; % For fast computation of correlations without deboor interpolation
if isfield(par, 'obsts')
    sigma = [];
    obbounds = zeros(2,length(par.obsts));
    for obst = 1:length(par.obsts)
        obbounds(1,obst) = length(sigma)+1;
        sigma = [sigma linspace(0,1,round(par.obsts(obst).N*fr) )];
        obbounds(2,obst) = length(sigma);
    end
    colLow = cell(length(par.obsts),1);
    rlow = par.r;
    for obst = 1:length(par.obsts)
        colLow{obst} = par.obsts(obst).colltau;
    end
else
    sigma = linspace(0,1,round(par.N*fr));
end
sr = length(sigma);
prevToc = toc;

if exist('A','var') && isfield(par, 'obsts') && strncmpi(method,'col',3) && ~strcmp(method,'colRint')
    R = zeros(par.N,sr);
    % We go over the columns of R so that we can reuse j1, j2 and wi
    for roi = 1:sr
        if (toc-prevToc > tim(1) )
            prevToc = toc;            
            display(['k=' num2str(par.k) ': ' num2str(roi/sr,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(sr-roi)/roi) ', est. end ' datestr(tim(2)+prevToc*sr/roi/24/3600)])
        end
        obsin = find((roi >= obbounds(1,:)) & (roi <= obbounds(2,:)));
        [~, cliloc] = min(abs(par.obsts(obsin).colltau-sigma(roi))); % Index i in oldColl (=col idx of old A2) closest to roit(roi)
        % Find range of indices in oldColl which correspond to nonzero A2(:,cliloc)
        nzIdxs = find(A(:,cliloc+par.r(1,obsin)-1))';
        if ~isempty(nzIdxs)
            T1 = sigma(roi)-Tcor;
            T2 = sigma(roi)+Tcor;
            [j1loc, j2loc, noSplit] = bounds2Ind(par.obsts(obsin).colltau,T1,T2);
            j1 = j1loc+par.r(1,obsin)-1;
            j2 = j2loc+par.r(1,obsin)-1;
            if noSplit
                wi = chi(par.obsts(obsin).colltau(j1loc:j2loc),T1, T1+pd*Tcor, T2-pd*Tcor, T2,0);
                for i = nzIdxs
                    R(i,roi) = (wi.*A(i,j1:j2))*c(j1:j2);
                end
            else
                wi1 = chi(par.obsts(obsin).colltau(j1loc:par.obsts(obsin).N),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
                wi2 = chi(par.obsts(obsin).colltau(1:j2loc),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1);
                for i = nzIdxs
                    R(i,roi) = (wi1.*A(i,j1:par.r(2,obsin)))*c(j1:par.r(2,obsin)) + (wi2.*A(i,par.r(1,obsin):j2))*c(par.r(1,obsin):j2);
                end
             end
        end
    end
    if strcmp(method,'colB')
        allBounds = zeros(4,1,par.N);
        % New columns will be added if more than one window is needed
        for i = 1:par.N
            obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
            tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
            [~, cli] = min(abs(colLow{obsin}-tc));
            cli = cli + rlow(1,obsin)-1;
            curThr = par.xi*max(abs(R(cli,:)));
            
            for obst = 1:length(par.obsts)
                rowt = R(cli,:);
                rowt(1:(obbounds(1,obst)-1)) = 0;
                rowt((obbounds(2,obst)+1):end) = 0;
                I = find(abs(rowt) >= curThr);
                ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
                collx = [];
                if (obst == obsin) && (numel(I) > 1)
                    bounds = [tc-a2w(1)/par.k, sigma(I(ini(1:(length(ini)-1) )+1)); tc+a2w(1)/par.k, sigma(I(ini(2:length(ini))))];
                elseif (numel(I) > 1)
                    bounds = [sigma(I(ini(1:(length(ini)-1) )+1)); sigma(I(ini(2:length(ini))))];
                    collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
                end
                decay = repmat(a2w(2)*a2w(1)/par.k,1,size(bounds,2));
                bounds = sortBounds(i,par,[(bounds + [-decay; decay]); [1;1]*decay], collx);
                allBounds(:,1:size(bounds,2),i) = bounds;
            end
        end
        % When making A2, do
%         for i = 1:par.N
%             if (toc-prevToc > printtoc)
%                 prevToc = toc;
%                 display(['Obst ' num2str(obstacle) ' at k=' num2str(ks(ki)) ': ', num2str(i/par.N,'%7.3f') '=A2%, now=' datestr(now)...
%                     ', est. # sec. left for A2=' num2str(toc*(par.N-i)/i) ', tmp comprErr=' num2str(norm(A2(1:(i-1),:)*c1-b(1:(i-1)))/...
%                     norm(b(1:(i-1))  )) ', est. end ' datestr(expectedEnd-(sum(v.timeA(:,2))/extf*(ppw(oi)*ks(ki)).^powTime - ...
%                     toc*par.N/i )/24/3600)]);
%             end
%             obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
%             collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
%             for obst = 1:length(par.obsts)
%                 if obst == obsin
%                     A2(i,par.r(1,obsin):par.r(2,obsin)) = windRow(i-par.r(1,obsin)+1,par.obsts(obsin),bounds,0,1);
%                 else
%                     A2(i,par.r(1,obst):par.r(2,obst)) = windRow('error',par.obsts(obst),bounds,0,1, collx);
%                 end
%             end
%         end        
        R = allBounds;
    end
    return
elseif exist('A','var') && ~isfield(par, 'obsts') && strncmpi(method,'col',3) && ~strcmp(method,'colRint')
    R = zeros(par.N, sr);
    % We go over the columns of R so that we can reuse j1, j2 and wi
    for roi = 1:sr
        if (toc-prevToc > tim(1) )
            prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(roi/sr,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(sr-roi)/roi) ', est. end ' datestr(tim(2)+prevToc*sr/roi/24/3600)])
        end
        [~, cliroi] = min(abs(par.colltau-sigma(roi))); % Index i in oldColl (=col idx of old A2) closest to roit(roi)
        % Find range of indices in oldColl which correspond to A2(:,cli) nonzero
        oldNzIdxs = find(A(:,cliroi))';
        if ~isempty(oldNzIdxs)
            T1 = sigma(roi)-Tcor;
            T2 = sigma(roi)+Tcor;
            [j1, j2, noSplit] = bounds2Ind(par.colltau,T1,T2);
            if noSplit
                wic = transpose(chi(par.colltau(j1:j2),T1, T1+pd*Tcor, T2-pd*Tcor, T2,0)).*c(j1:j2);
                for i = oldNzIdxs
                    R(i,roi) = A(i,j1:j2)*wic;
                end
            else
                wic1 = transpose(chi(par.colltau(j1:par.N),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1)).*c(j1:par.N);
                wic2 = transpose(chi(par.colltau(1:j2),T1, T1+pd*Tcor, T2-pd*Tcor, T2,1)).*c(1:j2);
                for i = oldNzIdxs
                    R(i,roi) = A(i,j1:par.N)*wic1 + A(i,1:j2)*wic2;
                end
            end
        end
    end
    
    if strcmp(method, 'colB')
        allBounds = zeros(4,1,par.N);
        curThr = par.xi*max(max(abs(R))); % Global threshold
        for i = 1:par.N
            if (toc-prevToc > tim(1))
                prevToc = toc;
                display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f') '=allBounds%, now=' datestr(now) ...
                    ', est. # sec. left for allBounds=' num2str(toc*(par.N-i)/i) ...
                    ', est. end ' datestr(tim(2) + prevToc*i/par.N/24/3600)]);
            end
            tc = par.colltau(i);
            [~, cli] = min(abs(par.colltau-tc)); % Index i in R closest to tc
            curThr = par.xi*max(abs(R(cli,:))); % Change to local threshold
            I = find(abs(R(cli,:)) >= curThr);
            ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
            bounds = [tc-a2w(1)/par.k, sigma(I(ini(1:(length(ini)-1) )+1)); tc+a2w(1)/par.k, sigma(I(ini(2:length(ini))))]; % id 1 on tc\pm T
            decay = repmat(a2w(2)*a2w(1)/par.k,1,size(bounds,2));
            bounds = sortBounds(i,par,[(bounds + [-decay; decay]); [1;1]*decay],[]); % size of bounds can change in this
            allBounds(:,1:size(bounds,2),i) = bounds;
        end
        R = allBounds;
    end
    return
elseif exist('A','var') && ~isfield(par, 'obsts') && strcmp(method,'rowB')
    allBounds = zeros(4,1,par.N);
    u = [par.t(end-par.dbf:end-1)-1, par.t, par.t(2:par.dbf+1)+1]; % the extended knots
    if fr == 1
        c1ip = [c; c(1:par.dbf)];
    else
        c1ip = deboor(par.dbf, u, [c; c(1:par.dbf)], sigma);
%         warning('Change code to use deboor interpolation for A');
    end
    Rrow = 0*sigma;
    jNs = nan*zeros(3,sr);
    for roi = 1:sr
       [j1, j2, noSplit] = bounds2Ind(sigma,sigma(roi)-Tcor,sigma(roi)+Tcor);
       jNs(:,roi) = [j1; j2; noSplit];
    end
    roi = round(sr/2); 
    [j1, j2] = bounds2Ind(sigma,sigma(roi)-Tcor,sigma(roi)+Tcor); % noSplit will be true
    wi = chi(sigma(j1:j2),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,0);
    
    for i = 1:par.N
        if (toc-prevToc > tim(1))
            prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f') '=allBounds%, now=' datestr(now) ...
                ', est. # sec. left for allBounds=' num2str(toc*(par.N-i)/i) ...
                ', est. end ' datestr(tim(2) + prevToc*i/par.N/24/3600)]);
        end
        tc = par.colltau(i);
        
        if fr == 1
            integrandt = [transpose(A(i,:)); transpose(A(i,1:par.dbf))].*c1ip;
        else
            integrandt = transpose(deboor(par.dbf,u, [transpose(A(i,:)); transpose(A(i,1:par.dbf))], sigma).*c1ip); % Use this if fr ~= 1
            % The variant above and below agree up to two digits and they are about as fast: fastest would be fr = 1
            % integrand = transpose(deboor(par.dbf,u, [transpose(A(i,:)).*c; transpose(A(i,1:par.dbf)).*c(1:par.dbf)], sigma));
        end
        
        nzCols = find(A(i,:)); % Will not be empty
        Rrow(:) = 0; % Need to reinitialise for each i
        for roi = nzCols
            if jNs(3,roi)
                Rrow(roi) = wi*integrandt(jNs(1,roi):jNs(2,roi)); % Faster than sum(wi.*integrand) and reuse wi
            else
                wi1 = chi(sigma(jNs(1,roi):sr), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                wi2 = chi(sigma(1:jNs(2,roi)), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                Rrow(roi) = wi1*integrandt(jNs(1,roi):sr) + wi2*integrandt(1:jNs(2,roi));
            end
        end
        
        curThr = par.xi*max(abs(Rrow)); % Local threshold
        I = find(abs(Rrow) >= curThr);
        ini = [0 find(I(2:end) - I(1:end-1) ~= 1) length(I)];
        bounds = [(tc-0.02), sigma(I(ini(1:(length(ini)-1) )+1)); (tc+0.02), sigma(I(ini(2:length(ini))))];
        % Or: bounds = [tc-a2w(1)/par.k, sigma(I(ini(1:(length(ini)-1) )+1)); tc+a2w(1)/par.k, sigma(I(ini(2:length(ini))))]; % id 1 on tc\pm T
        decay = repmat(0.02,1,size(bounds,2));
        decay(1) = 0.02*16/par.k;
        % Or: decay = repmat(a2w(2)*a2w(1)/par.k,1,size(bounds,2));
        bounds = sortBounds(i,par,[(bounds + [-decay; decay]); [1;1]*decay], []); % size of bounds can change in this
        allBounds(:,1:size(bounds,2),i) = bounds;
    end
    R = allBounds;
    return
end

% Compute the correlations as an integral through a simple Riemann-sum without normalisation.
% We could compute R column by column because then the window around sigma(roi) would be constant, so that we do not have to
% recompute the same window indices j1 and j2 for each row. But this would need to recompute the integrand with the
% (expensive) Bessel function for each column.
if ~isfield(par, 'obsts') % Single scattering obstacle without using A
    R = zeros(par.N, sr);
    u = [par.t(end-par.dbf:end-1)-1, par.t, par.t(2:par.dbf+1)+1]; % the extended knots
    c1ip = deboor(par.dbf, u, [c; c(1:par.dbf)], sigma);
    
    % Find the intervals of u where the respective sigma lies.
    for i = 1:par.N
        if (toc-prevToc > tim(1) )
            prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(par.N-i)/i) ', est. end ' datestr(tim(2)+prevToc*par.N/i/24/3600)])
        end
        tc = par.colltau(i);
        integrand = 1i/4.*besselh(0, 1, par.k*sqrt(sum((par.par(tc*ones(size(sigma))) ...
            -par.par(sigma) ).^2, 1)) ).*c1ip.*par.gradnorm(sigma);
        integrand(isnan(integrand)) = 0;
        for roi = 1:size(R,2)
            [j1, j2, noSplit] = bounds2Ind(sigma,sigma(roi)-Tcor,sigma(roi)+Tcor);
            if noSplit
                wi = chi(sigma(j1:j2),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,0);
                R(i,roi) = sum(wi.*integrand(j1:j2) );
            else
                wi1 = chi(sigma(j1:sr),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                wi2 = chi(sigma(1:j2), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                R(i,roi) = sum(wi1.*integrand(j1:sr) ) + sum(wi2.*integrand(1:j2) );
            end
        end
    end
    return
end
% Multiple scattering obstacle without using A
R = zeros(par.N,sr);
c1ip = NaN*sigma;
for obst = 1:length(par.obsts)
    u = [par.obsts(obst).t(end-par.obsts(obst).dbf:end-1)-1, par.obsts(obst).t, par.obsts(obst).t(2:par.obsts(obst).dbf+1)+1]; % the extended knots
    c1ip(obbounds(1,obst):obbounds(2,obst)) = deboor(par.obsts(obst).dbf, u, ...
        [c(par.r(1,obst):par.r(2,obst)); c(par.r(1,obst):(par.r(1,obst) +par.obsts(obst).dbf-1))], sigma(obbounds(1,obst):obbounds(2,obst)) );
end
for i = 1:par.N
    if (toc-prevToc > tim(1) )
        prevToc = toc;
            display(['k=' num2str(par.k) ': ' num2str(i/par.N,'%7.3f')  '=R%, now=' datestr(now) ', est. # sec. left for R=' ...
                num2str(toc*(par.N-i)/i) ', est. end ' datestr(tim(2)+prevToc*par.N/i/24/3600)])
    end
    obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
    tc = par.obsts(obsin).colltau(i-par.r(1,obsin)+1);
    colx = par.obsts(obsin).par(tc);
    for obst = 1:length(par.obsts)
        cusigma = sigma(obbounds(1,obst):obbounds(2,obst));
        integrand = 1i/4.*besselh(0, 1, par.k*sqrt(sum((repmat(colx,1,length(cusigma)) -par.obsts(obst).par(cusigma) ).^2, 1)) ).*...
            c1ip(obbounds(1,obst):obbounds(2,obst)).*par.obsts(obsin).gradnorm(cusigma);
        integrand(isnan(integrand)) = 0;
        for roi = obbounds(1,obst):obbounds(2,obst)
            [j1, j2, noSplit] = bounds2Ind(cusigma,sigma(roi)-Tcor,sigma(roi)+Tcor);
            if noSplit
                wi = chi(cusigma(j1:j2),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,0);
                R(i,roi) = sum(wi.*integrand(j1:j2) );
            else
                wi1 = chi(cusigma(j1:length(cusigma)),sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                wi2 = chi(cusigma(1:j2), sigma(roi)-Tcor, sigma(roi)-(1-pd)*Tcor, sigma(roi)+Tcor*(1-pd),sigma(roi)+Tcor,1);
                R(i,roi) = sum(wi1.*integrand(j1:length(cusigma)) ) + sum(wi2.*integrand(1:j2) );
            end
        end
    end
    debug = 1;
end

end

% Evaluate the spline of degree k determined by the knots t(n+2*k+1) and coefficients c(n+k) at the points in x.
function s = deboor(k, t, c, x)

s = zeros(size(x));
for in=1:length(x)
    % Determine the interval where x(in) lies
    if x(in) == t(length(t)-k)
        j = length(t) - 2*k - 2;
    else
        j = find(x(in) >= t, 1, 'last') - k - 1;
    end
    d=c(j-k+k+1:j+k+1);
    for r=1:k
        for i=j:-1:j-k+r
            idx = i+k+1;
            alfa = (x(in) - t(idx)) / (t(idx+k+1-r)-t(idx));
            d(idx-j) = alfa*d(idx-j) + (1-alfa)*d(idx-j-1);
        end
    end
    s(in) = d(k+1);
end

end
