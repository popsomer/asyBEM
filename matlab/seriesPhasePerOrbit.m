% Compute the series expansion of the phase of the periodic orbit mode
% Input
%	par       - The par structure, containing an 'obsts' field in which the J obstacles are ordered according to the periodic orbit. When J > 2,
%       this function should be called with all collections of J-1 obstacles, J-2 and so on till J-J+2 to obtain the phases of all
%       periodic orbits. ?? The order is important, but supplying par.obsts(J:-1:1) should give the same results ?? 
%   maxOrder  - The maximal order of the series expansion
% Output
%   taus  - The J parameters on the respective obstacles giving the periodic orbit
%   c     - The expansion coefficients of the phase
%   a     - The expansion coefficients of the stationary point
%   ft    - Series expansion coefficients of the distance
function [taus, c, a, ft] = seriesPhasePerOrbit(par, maxOrder)
% Test using par = getObst(5); maxOrder = 5; [taus, c, a, ft] = seriesPhasePerOrbit(par, maxOrder)
J = length(par.obsts);

% Compute a good initial guesses for the taus
midpts = zeros(2,J);
taut = linspace(0, 1, 100); 
taut(end) = [];
for obst = 1:J
    midpts(:,obst) = mean(par.obsts(obst).par(taut), 2);
end
mid = repmat(mean(midpts, 2), 1, length(taut));

taus = nan(J,1);
for obst = 1:J
    partmp = par.obsts(obst).par(taut);
    [~, ix] = min( (partmp(1,:) -mid(1,:) ).^2 + (partmp(2,:) -mid(2,:) ).^2);
    taus(obst) = taut(ix);
end

% Improve the initial guesses
    function dt = distmp(ts)
        dt = norm(par.obsts(1).par(ts(1)) -par.obsts(end).par(ts(end)) );
        for ob = 2:length(par.obsts)
            dt = dt + norm(par.obsts(ob-1).par(ts(ob-1)) -par.obsts(ob).par(ts(ob)) );
        end
    end
taus = fminsearch(@distmp, taus);

binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n)*0^((n<0) + abs(rem(n,1)) ); % Not to be used vectorized (and n should be a positive integer)
% ft = nan(J, maxOrder, maxOrder);
ft = nan(J, maxOrder+1, maxOrder+1);% Actually, we only need (1, maxOrder, maxOrder) if J = 2
d = nan(J,1);
for obst = 1:J
    other = obst-1;
    if other == 0
        other = J;
    end
%     Gama = par.obsts(obst).serpar(taus(obst), maxOrder);
    Gama = par.obsts(obst).serpar(taus(obst), maxOrder+1);
    Gamb = par.obsts(other).serpar(taus(other), maxOrder+1);
    Lambda = nan(maxOrder+1, maxOrder+1, 2);
    Lambda(1,1,:) = Gamb(:,1).^2 - 2*Gama(:,1).*Gamb(:,1) + Gama(:,1).^2;
    d(obst) = sqrt(sum(Lambda(1,1,:)));
    for l = 2:maxOrder+1
        Lambda(1,l,:) = sum(Gamb(:,1:l).*Gamb(:, l:-1:1), 2) -2*Gama(:,1).*Gamb(:,l);
    end
    for j = 2:maxOrder+1
        Lambda(j,1,:) = sum(Gama(:,1:j).*Gama(:, j:-1:1), 2) -2*Gama(:,j).*Gamb(:,1);
    end
    for l = 2:maxOrder+1
        for j = 2:maxOrder+1
            Lambda(j,l,:) = 2*Gama(:,j).*Gamb(:,l);
        end
    end
%     zt = nan(2*maxOrder-2, maxOrder, maxOrder);
    zt = nan(2*maxOrder, maxOrder+1, maxOrder+1);
    zt(1,:,:) = (Lambda(:,:,1) + Lambda(:,:,2))/(Lambda(1,1,1) + Lambda(1,1,2));
    zt(1,1,1) = nan;
    for m = 2:2*maxOrder
        for k = 1:maxOrder+1
            for i = max(1, m-k+2):maxOrder+1
                zt(m,i,k) = 0;
                for r = 1:k
                    s = (1 + (r==1)):(i-max(0,m-k+r-1));
                    zt(m,i,k) = zt(m,i,k) + sum(zt(1,s,r).*zt(m-1, i-s+1, k-r+1)); 
                end
            end
        end
    end
%     for k = 1:maxOrder
    for k = 1:maxOrder+1
        for i = (1+(k==1)):maxOrder+1
            ft(obst, i,k) = zt(1,i,k)/2;
            for m = 2:k-2+i
                ft(obst,i,k) = ft(obst,i,k) + binom(1/2,m)*zt(m,i,k);
            end
            ft(obst,i,k) = d(obst)*ft(obst,i,k);
        end
    end
end

a = nan(maxOrder, maxOrder, J);
c = nan(maxOrder, J);
da = 0; % Shift in a due to undefined a(1,1)
if J == 2
    da = 1;
%     if (ft(1, 2, 1) ~= 0) || (ft(1, 1, 2) ~= 0)
    if (abs(ft(1, 2, 1)) > d(1)*20*eps) || (abs(ft(1, 1, 2)) > d(1)*20*eps)
        error('Distance should be minimal.');
    elseif (abs(ft(1, 3, 1)) < d(1)*20*eps) && (abs(ft(1, 1, 3)) < d(1)*20*eps) && (abs(ft(1, 2,2)) < d(1)*20*eps)
        da = 2;
        warning('Cubic equation is not implemented.');
    elseif (abs(ft(1, 3, 1)) < d(1)*20*eps) || (abs(ft(1, 1, 3)) < d(1)*20*eps) || (abs(ft(1, 2,2)) < d(1)*20*eps)
        error('No limit cycle exists.');
    elseif (ft(1,3,1) < 0) || (ft(1,1,3) < 0)
        error('Distance should be minimal.')
    end
    a(1,1,1) = -2*ft(1,1,3)/ft(1,2,2) + sign(ft(1,2,2))*sqrt((2*ft(1,1,3)/ft(1,2,2))^2 - ft(1,1,3)/ft(1,3,1));
%     a(1,1,2) = -2*ft(1,3,1)/ft(1,2,2) + sign(ft(1,2,2))*sqrt((2*ft(1,3,1)/ft(1,2,2))^2 - ft(1,3,1)/ft(1,1,3));
    a(1,1,2) = -ft(1,2,2)/(a(1,1,1)*ft(1,1,2) +4*ft(1,1,3));
    c(1,1) = 0;
    c(1,2) = 0;
    c(2,1) = -ft(1,2,2)/2/a(1,1,1) - ft(1,3,1);
    c(2,2) = -ft(1,2,2)/2/a(1,1,2) - ft(1,1,3);
else
    c(1,1) = -ft(1,2,1)/2;
%     a(1,1,1) = 2*(-ft(J,2,1)/2)/ft(1,2,1) -2;
    a(1,1,1) = -ft(J,2,1)/ft(1,2,1) -2;
    for obst = 2:J
        c(1,obst) = -ft(obst,2,1)/2;
        a(1,1,obst) = 2*c(1,obst-1)/ft(obst,2,1) -2;
    end
end

nc = norm(c(1:1+da,:) );
na = norm(squeeze(a(1,1,:)) );
if (nc == inf) || isnan(nc) || (nc == 0) || (na == inf) || isnan(na) || (na == 0)
    error('Incorrect results');
end

for i=2:maxOrder
    for obst = 1:J
        for ii = 2:i
            for l = ii:maxOrder
                k = (ii-1):(l-1);
                a(ii,l,obst) = sum(a(ii-1,k,obst).*a(1,l-k,obst) );
            end
        end
        c(i+da, obst) = 0;
        for l = 1:i-1
            n = 3:(i-l+2);
%             c(i+da, obst) = c(i+da, obst) + sum((n-1).*ft(obst,n,l).*a(n-2,i-l,obst));
            c(i+da, obst) = c(i+da, obst) + (n-1).*ft(obst,n,l)*a(n-2,i-l,obst);
        end
        k = 2:i-1;
        c(i+da, obst) = -1/i/a(i-1,i-1,obst)*(ft(obst,2,i) + c(i+da, obst) + k*(c(k,obst).*a(k-1,i-1,obst)) );
    end
    
    for obst = 1:J
        other = obst-1;
        if other == 0
            other = J;
        end
        a(1,i,obst) = 0;
        tmp = 0;
        for l = 1:i-1
            n = 3:(i-l+2);
            a(1,i,obst) = a(1,i,obst) + (ft(obst,n,l).*(n-1))*a(n-2,i-l,obst);
            n = l+1;
            m = (n-1):(i-1);
            tmp = tmp + a(n-1,m)*squeeze(ft(obst,n,i+1-m));
        end
        k = 2:(i-1);
        r = 3:(i+1);
        a(1,i,obst) = (a(1,1,obst)/i*(ft(obst,2,i) + a(1,i,obst) + sum(c(k,obst).*k'.*a(k-1,i-1,obst) ) ) ...
            -sum(c(k,obst).*a(k,i,obst)) + c(i,other) -ft(obst,r,1)*a(r-1,i,obst) -ft(obst,1,i+1) -tmp)/(c(1,obst) + ft(obst,2,1));
    end
    
%     na = norm(a(1:i,:) );
    nc = norm(c(1:1+da,:) );
    na = norm(squeeze(a(1,1:i,:)) );
    if (nc == inf) || isnan(nc) || (nc == 0) || (na == inf) || isnan(na) || (na == 0)
        error(['Incorrect results at i = ' num2str(i)]);
    end
end
    
end
