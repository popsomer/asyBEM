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
%   d     - Distances between the obstacles along the periodic orbit
function [taus, c, a, ft, d] = seriesPhasePerOrbit(par, maxOrder)
% Test using par = getObst(5); maxOrder = 5; [taus, c, a, ft] = seriesPhasePerOrbit(par, maxOrder)
% [norm(c(1,:)), norm(c(2,:)-sqrt(2)*pi^2), norm(c(3,:)), norm(c(4,:) +11/12*sqrt(2)*pi^4), norm(c(5,:)), norm(c(6,:) - 2783/2520*sqrt(2)*pi^6), norm(c(7,:)), norm(c(8,:) + 358021/205632*sqrt(2)*pi^8)]
% [norm(squeeze(a(1,1,:) - (3-2*sqrt(2)))), norm(squeeze(a(1,2,:))), norm(squeeze(a(1,3,:)-7*(17*sqrt(2)-24)*pi^2)), norm(squeeze(a(1,4,:))), norm(squeeze(a(1,5,:) -(1205811*sqrt(2) -1705312)*pi^4/84)), norm(squeeze(a(1,6,:))), norm(squeeze(a(1,7,:)  + pi^6/128520*(289615597399*sqrt(2) -409578202752)))]

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
% taus = fminsearch(@distmp, taus);
% accTaus = fminsearch(@(ts) norm(par.obsts(1).par(ts(1)) -par.obsts(2).par(ts(2)) ), taus, optimset('Display', 'iter', 'TolX', sqrt(eps)) )
% taus = fminsearch(@distmp, taus, optimset('TolX', eps) );

function [dt, g] = distWgrad(ts)
    dt = 0;
    g = zeros(length(ts), 1);
    prev = length(ts);
    prGam = par.obsts(prev).serpar(ts(prev),2);
%     prGam = par.obsts(end).serpar(ts(end),1);
    for ob = 1:length(par.obsts)
        gam = par.obsts(ob).serpar(ts(ob), 2);
%         difGam = gam -prGam;
        difGam = gam(:,1) -prGam(:,1);
        curDis = norm(difGam);
        dt = dt + curDis;
%         g(ob,:) = transpose(difGam(:,1))*gam(:,2)/curDis;
        g(ob) = g(ob) + transpose(difGam)*gam(:,2)/curDis;
        g(prev) = g(prev) - transpose(difGam)*prGam(:,2)/curDis;
        prev = ob;
        prGam = gam;
    end
end
% taus = fminunc(@distWgrad, taus, optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true, 'TolX', eps));
taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, ...
    'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000));

% Try to improve to machine precision using a few Newton iterations?? or Hessian:

function [dt, g, H] = distWHess(ts)
    dt = 0;
    g = zeros(length(ts), 1);
    H = zeros(length(ts), length(ts));
    prev = length(ts);
    prGam = par.obsts(prev).serpar(ts(prev), 3);
    for ob = 1:length(par.obsts)
        gam = par.obsts(ob).serpar(ts(ob), 3);
        difGam = gam(:,1) -prGam(:,1);
        curDis = norm(difGam);
        dt = dt + curDis;
        gob = transpose(difGam)*gam(:,2)/curDis;
        g(ob) = g(ob) + gob;
        gprev = - transpose(difGam)*prGam(:,2)/curDis;
        g(prev) = g(prev) + gprev;
        
        H(ob,ob) = H(ob,ob) +(sum(gam(:,3).^2) +transpose(difGam)*gam(:,2).^2 -g(ob)/2/curDis)/curDis;
        H(ob,prev) = H(ob,prev) +(sum(gam(:,2).*prGam(:,2)) - gprev*gob*0.9 )/curDis;
%         H(ob,prev) = H(ob,prev) +(sum(gam(:,2).*prGam(:,2)) - gprev*gob )/curDis;
        H(prev,ob) = H(prev,ob) +(sum(gam(:,2).*prGam(:,2)) - gob*gprev )/curDis;
%         H(prev,ob) = H(prev,ob) +(sum(gam(:,2).*prGam(:,2)) + g(ob)*g(prev) )/curDis;
        H(prev,prev) = H(prev,prev)  -(sum(prGam(:,3).^2) +transpose(difGam)*prGam(:,2).^2*1.1 + g(prev)/2/curDis)/curDis;
%         H(prev,prev) = H(prev,prev)  -(sum(prGam(:,3).^2) +transpose(difGam)*prGam(:,2).^2 + g(prev)/2/curDis)/curDis;
%         if abs(H(ob,prev) - H(prev,ob) ) > 20*eps
%             error('Hessian not symmetric');
%         end
        prev = ob;
        prGam = gam;
    end
%     if norm(transpose(H) -H) > 20*eps*norm(H)
%         error('Hessian not symmetric')
%     end
end

% taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true, ...
% taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'newton', 'SpecifyObjectiveGradient', true, ...
% taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, ...
%     'HessianFcn','objective', 'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000));
%'Algorithm', 'trust-region' % 'interior-point'
% [dt, g, H] = distWHess(taus) % check optimality

binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n)*0^((n<0) + abs(rem(n,1)) ); % Not to be used vectorized (and n should be a positive integer)
% ft = nan(J, maxOrder, maxOrder);
if 0
    fj = J - (J==2);
else
    fj = J;
end
% ft = nan(J, maxOrder+1, maxOrder+1);
% Actually, we only need (1, maxOrder, maxOrder) if J = 2 but this would slightly complicate the implementation.
ft = nan(fj, maxOrder+1, maxOrder+1);
d = nan(fj,1);
% d = nan(J,1);
for obst = 1:fj
% for obst = 1:J
%     other = obst-1;
%     if other == 0
%         other = J;
%     end
    % Above takes prev, below takes next as in makeV1.m
    other = obst +1 -(obst == length(par.obsts))*length(par.obsts);
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
%             Lambda(j,l,:) = 2*Gama(:,j).*Gamb(:,l);
            Lambda(j,l,:) = -2*Gama(:,j).*Gamb(:,l);
        end
    end
%     zt = nan(2*maxOrder-2, maxOrder, maxOrder);
    zt = nan(2*maxOrder, maxOrder+1, maxOrder+1);
%     zt(1,:,:) = (Lambda(:,:,1) + Lambda(:,:,2))/(Lambda(1,1,1) + Lambda(1,1,2));
    % Avoid wrong conversion of dimensions by this loop
    for i = 1:maxOrder+1
        zt(1,i,:) = (Lambda(i,:,1) + Lambda(i,:,2))/d(obst)^2;
%         for k = 1:maxOrder+1
%         end
    end
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

a = nan(maxOrder, maxOrder, J); % a(1,maxOrder,:) should not be used when J==2
% a = nan(maxOrder -(J==2), maxOrder -(J==2), J);
c = nan(maxOrder, J);
% da = 0; % Shift in a due to undefined a(1,1)
if J == 2
%     da = 1;
%     if (ft(1, 2, 1) ~= 0) || (ft(1, 1, 2) ~= 0)
%     if (abs(ft(1, 2, 1)) > d(1)*20*eps) || (abs(ft(1, 1, 2)) > d(1)*20*eps)
%     tol = d(1)*1e-10; 
    tol = d(1)*1e-7;
    % No tolerance 20*eps but 1e-10 due to loss of digits in objective fct in fmincon for quasi-Newton
    % with exact parametrisation without FFT. Even 1e-6*d(1) for trust-region which uses the analytic Hessian.
    if (abs(ft(1, 2, 1)) > tol) || (abs(ft(1, 1, 2)) > tol) 
        error('Distance should be minimal.');
%     elseif (abs(ft(1, 3, 1)) < d(1)*20*eps) && (abs(ft(1, 1, 3)) < d(1)*20*eps) && (abs(ft(1, 2,2)) < d(1)*20*eps)
    elseif (abs(ft(1, 3, 1)) < tol) && (abs(ft(1, 1, 3)) < tol) && (abs(ft(1, 2,2)) < tol)
%         da = 2;
        warning('Cubic equation is not implemented.');
%     elseif (abs(ft(1, 3, 1)) < d(1)*20*eps) || (abs(ft(1, 1, 3)) < d(1)*20*eps) || (abs(ft(1, 2,2)) < d(1)*20*eps)
    elseif (abs(ft(1, 3, 1)) < tol) || (abs(ft(1, 1, 3)) < tol) || (abs(ft(1, 2,2)) < tol)
        error('No limit cycle exists.');
    elseif (ft(1,3,1) < 0) || (ft(1,1,3) < 0)
        error('Distance should be minimal.')
    end
    a(1,1,1) = -2*ft(1,1,3)/ft(1,2,2) + sign(ft(1,2,2))*sqrt((2*ft(1,1,3)/ft(1,2,2))^2 - ft(1,1,3)/ft(1,3,1));
%     a(1,1,2) = -2*ft(1,3,1)/ft(1,2,2) + sign(ft(1,2,2))*sqrt((2*ft(1,3,1)/ft(1,2,2))^2 - ft(1,3,1)/ft(1,1,3));
%     a(1,1,2) = -ft(1,2,2)/(a(1,1,1)*ft(1,1,2) +4*ft(1,1,3));
    a(1,1,2) = -ft(1,2,2)/(a(1,1,1)*ft(1,2,2) +4*ft(1,1,3));
%     for i = 2:maxOrder-1 %maxOrder
%         a(i,i,:) = a(i-1,i-1,:).*a(1,1,:);
%     end
    c(1,1) = 0;
    c(1,2) = 0;
    c(2,1) = -ft(1,2,2)/2/a(1,1,1) - ft(1,3,1);
    c(2,2) = -ft(1,2,2)/2/a(1,1,2) - ft(1,1,3);
    for i = 3:maxOrder
        for l = (i-1):maxOrder % :maxOrder-1
            k = (l-i+2):(l-1);
            a(l-i+3,l,:) = sum(a(l-i+2,k,:).*a(1,l-k,:), 2);
        end
        M = zeros(4,4);
        M(1,1) = a(i,i,1);
        M(1,2) = -1; 
%         M(1,3) = ft(1,2,2) +2*ft(1,3,1)*a(1,1,1) + 2*c(2,1);
%         M(1,3) = ft(1,2,2) +2*a(1,1,1)*(ft(1,3,1) + c(2,1));
        % Both M(1,3) and M(2,4) should be zero by construction
        M(2,2) = a(i,i,2);
        M(2,1) = -1;
%         M(2,4) = ft(1,2,2) +2*a(1,1,2)*(ft(1,1,3) + c(2,2));
        M(3,1) = i*a(i-1,i-1,1);
        M(3,3) = 2*(ft(1,3,1) + c(2,1));
        M(4,2) = i*a(i-1,i-1,2);
        M(4,4) = 2*(ft(1,1,3) + c(2,2));
        rhs = -[(ft(1,1,i+1)+ sum(c(3:i-1,1).*a(3:i-1,i,1))); ...
            (ft(1,i+1,1) + sum(c(3:i-1,2).*a(3:i-1,i,2))); ...
            (ft(1,2,i) +(3:i-1)*(c(3:i-1,1).*a(2:i-2,i-1,1))) ; ...%(ft(1,2,i) +(c(3:i-1,1).*(2:i-1) )*a(2:i-2,i-1,1)) ; ...
            (ft(1,i,2) +(3:i-1)*(c(3:i-1,2).*a(2:i-2,i-1,2)))];
        rhs(1) = rhs(1) - (ft(1,3,1) + c(2,1))*sum(a(1,2:i-2,1).*a(1,i-2:-1:2,1)) - ft(1,4:i+1,1)*a(3:i,i,1);
        rhs(2) = rhs(2) - (ft(1,1,3) + c(2,2))*sum(a(1,2:i-2,2).*a(1,i-2:-1:2,2)) - sum(squeeze(ft(1,1,4:i+1)).*a(3:i,i,2));
        for l=1:i-1
            rhs(1) = rhs(1) - a(l,l:(i-1-(l==1)),1)*squeeze(ft(1,l+1,(i+1-l):-1:(2+(l==1)) ) ); % i+1-(l:i-1-(l==1))));
            rhs(2) = rhs(2) - sum(a(l,l:(i-1-(l==1)),2).*ft(1,(i+1-l):-1:(2+(l==1)), l+1) );
%             rhs(3) = rhs(3) -(ft(1,(3+(l==1)):(i-l+2),l).*((3+(l==1)):(i-l+2)))*a( (1+(l==1)):(i-l), i-l,1);
%             rhs(4) = rhs(4) -((3+(l==1)):(i-l+2))*(squeeze(ft(1,l,(3+(l==1)):(i-l+2))).*a( (1+(l==1)):(i-l), i-l,2));
            rhs(3) = rhs(3) -(ft(1,(3+(l==1)):(i-l+2),l).*((2+(l==1)):(i-l+1)))*a( (1+(l==1)):(i-l), i-l,1);
            rhs(4) = rhs(4) -((2+(l==1)):(i-l+1))*(squeeze(ft(1,l,(3+(l==1)):(i-l+2))).*a( (1+(l==1)):(i-l), i-l,2));
        end
%         s = rhs\M;
        s = M\rhs;
        c(i,:) = s(1:2);
        a(1,i-1,:) = s(3:4);
    end
%     return
% else
end

% c(1,1) = -ft(1,2,1)/2;
%     a(1,1,1) = 2*(-ft(J,2,1)/2)/ft(1,2,1) -2;
% a(1,1,1) = -ft(J,2,1)/ft(1,2,1) -2;
% for obst = 2:J
%     c(1,obst) = -ft(obst,2,1)/2;
%     a(1,1,obst) = 2*c(1,obst-1)/ft(obst,2,1) -2;
% end
for obst = 1:J
    c(1,obst) = -ft(obst,2,1);
end
if (J == 2) && (norm(c(1,:) ) > tol)
    error('taus does not give the minimum of the distance.')
elseif (J ~= 2) && (norm(ft(:,1,2) -transpose([c(1,2:end), c(1,1)])) > norm(c(1,:))*eps^(3/7))
    error('taus does not give the minimum of the distance.')
end


    % op = [omega^alpha, omega^beta, ... ; psi^alpha, ...]
    % ca = [c^alpha, c^beta, ... ; a^alpha, ... ]
    function op = omPsi(ca) 
        op = [[-ca(1,2:end), -ca(1,1)]; 2*ca(1,:).*ca(2,:)];
        for ob = 1:size(ca,2)
%             op(1,ob) = op(1,ob) +ft(ob,1,3) - ft(ob,2,2)/2*ca(2,ob);
            op(1,ob) = op(1,ob) +ft(ob,1,3) + ft(ob,2,2)/2*ca(2,ob);
            op(2,ob) = op(2,ob)  + ft(ob,2,2) +2*ft(ob,3,1)*ca(2,ob);
        end
    end

ca = fsolve(@omPsi, [+transpose(ft(:,3,1)); zeros(1,J)]);
if (norm(c(2,:)-ca(1,:)) > eps^(1/2)*norm(c(2,:)) ) || (norm(squeeze(a(1,1,:)) -transpose(ca(2,:))) > eps^(1/2)*norm(ca(2,:)) )
    error('Different result');
elseif any(c(2,:) < eps^(1/2))
    warning('c2 is zero or negative');
end
c(2,:) = ca(1,:);
a(1,1,:) = ca(2,:);

for i = 3:maxOrder
    for l = (i-1):maxOrder % :maxOrder-1
        k = (l-i+2):(l-1);
        a(l-i+3,l,:) = sum(a(l-i+2,k,:).*a(1,l-k,:), 2);
    end
    M = zeros(2*J, 2*J);
    rhs = zeros(2*J,1); %-[ft(:,1,i+1); ft(:,2,i)];
    for obst = 1:J
        other = obst +1 -(obst == J)*J;
        M(obst, obst) = a(i,i,obst);
        M(obst, other) = -1;
        % a_{i-1} should not influence omega_i because of omega_2
        M(J+obst, obst) = i*a(i-1,i-1,obst);
        M(J+obst, J+obst) = 2*(ft(obst,3,1) + c(2,obst));
        
        rhs(obst) = -sum(c(3:i-1,obst).*a(3:i-1,i,obst)) - ft(obst,4:i+1,1)*a(3:i,i,obst) ...
           -ft(obst,1,i+1) - (ft(obst,3,1) + c(2,obst))*sum(a(1,2:i-2,obst).*a(1,i-2:-1:2,obst));
        rhs(J+obst) = -ft(obst,2,i) - (3:i-1)*(c(3:i-1,obst).*a(2:i-2,i-1,obst));
        
        for l=1:i-1
            rhs(obst) = rhs(obst) - a(l, l:(i-1-(l==1)), obst)*squeeze(ft(obst,l+1,(i+1-l):-1:(2+(l==1)) ) );
            rhs(J+obst) = rhs(J+obst) -(ft(obst, (3+(l==1)):(i-l+2),l).*((2+(l==1)):(i-l+1)))*a( (1+(l==1)):(i-l), i-l, obst);
        end
    end
    s = M\rhs; 
    if (J == 2) && ( (norm(transpose(c(i,:)) -s(1:2)) > eps^(1/2)*norm(c(i,:)) ) || ...
            (norm(squeeze(a(1,i-1,:)) -s(3:4)) > eps^(1/2)*norm(s(3:4))) )
        error('Different result');
    end
    c(i,:) = s(1:J);
    a(1,i-1,:) = s(J+1:end);
end

return 
% Below old attempts
for i=2:maxOrder
    for obst = 1:J
%         for ii = 2:i
%             for l = ii:maxOrder
%                 k = (ii-1):(l-1);
%                 a(ii,l,obst) = sum(a(ii-1,k,obst).*a(1,l-k,obst) );
%             end
%         end
        for l = i:maxOrder
            k = (l-i+1):(l-1);
            a(l-i+2,l,obst) = sum(a(l-i+1,k,obst).*a(1,l-k,obst) );
        end
        c(i, obst) = 0;
        for l = 1:i-1
            n = 3:(i-l+2);
            c(i, obst) = c(i, obst) + (n-1).*ft(obst,n,l)*a(n-2,i-l,obst);
        end
        k = 2:i-1;
        c(i, obst) = -1/i/a(i-1,i-1,obst)*(ft(obst,2,i) + c(i, obst) + k*(c(k,obst).*a(k-1,i-1,obst)) );
    end
    
    for obst = 1:J
%         other = obst-1;
%         if other == 0
%             other = J;
%         end
        % Above takes previous, below takes next as in makeV1.m
        other = obst +1 -(obst == length(par.obsts))*length(par.obsts);
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
    
    nc = norm(c(1:i,:) );
    na = norm(squeeze(a(1,1:i,:)) );
    if (nc == inf) || isnan(nc) || (nc == 0) || (na == inf) || isnan(na) || (na == 0)
        error(['Incorrect results at i = ' num2str(i)]);
    end
end


% a = a(:,1:maxOrder-da,:); % Need a(k-1,k-1,obst) too

% nc = norm(c(1+da,:) );
% na = norm(squeeze(a(1,1,:)) );
% if (nc == inf) || isnan(nc) || (nc == 0) || (na == inf) || isnan(na) || (na == 0)
%     error('Incorrect results');
% end
% da = nan

% for i=2:maxOrder
%     
%  %         c(i, obst) = 0;
% %         for l = 1:i-1
% %             n = 3:(i-l+2);
% % %             c(i+da, obst) = c(i+da, obst) + sum((n-1).*ft(obst,n,l).*a(n-2,i-l,obst));
% %             c(i, obst) = c(i, obst) + (n-1).*ft(obst,n,l)*a(n-2,i-l,obst);
% %         end
% % %         k = 2:i-1;
% % %         c(i+da, obst) = -1/i/a(i-1,i-1,obst)*(ft(obst,2,i) + c(i+da, obst) + k*(c(k,obst).*a(k-1,i-1,obst)) );
% %         k = 2:i-1;
% %         c(i, obst) = -1/i/a(i-1,i-1,obst)*(ft(obst,2,i) + c(i, obst) + k*(c(k,obst).*a(k-1,i-1,obst)) );
%     
%     for obst = 1:J
%         other = obst-1;
%         if other == 0
%             other = J;
%         end
%         a(1,i,obst) = 0;
%         tmp = 0;
% %         for l = 1:i-1
%         for l = 1:i+da-1
% %             n = 3:(i-l+2);
%             n = (3 + da*(l==1)):(i+da-l+2);
% %             a(1,i,obst) = a(1,i,obst) + (ft(obst,n,l).*(n-1))*a(n-2,i-l,obst);
%             a(1,i,obst) = a(1,i,obst) + (ft(obst,n,l).*(n-1))*a(n-2,i+da-l,obst);
%             n = l+1;
% %             m = (n-1):(i-1);
%             m = (n-1):(i+da-1-da*(n==2));
% %             tmp = tmp + a(n-1,m)*squeeze(ft(obst,n,i+1-m));
%             tmp = tmp + a(n-1,m)*squeeze(ft(obst,n,i+da+1-m));
%         end
% %         k = 2:(i-1);
%         k = (2 + da):(i+da-1);
% %         r = 3:(i+1);
%         r = (3 + da):(i+da+1);
% %         a(1,i,obst) = (a(1,1,obst)/i*(ft(obst,2,i) + a(1,i,obst) + sum(c(k,obst).*k'.*a(k-1,i-1,obst) ) ) ...
% %             -sum(c(k,obst).*a(k,i,obst)) + c(i,other) -ft(obst,r,1)*a(r-1,i,obst) -ft(obst,1,i+1) -tmp)/(c(1,obst) + ft(obst,2,1));
%         denom = c(1,obst) + ft(obst,2,1) + da*( ft(obst,2,2) + 2*(c(2,obst) + ft(obst,3,1) )*(1 +a(1,1,obst)) );
%         extra = da*(c(2,obst) +ft(obst,3,1) )*sum(a(1,2:i+da-2,obst).*a(1,i+da-(2:i+da-2), obst));
%         if ~isempty(k)
%             extra = extra + a(1,1,obst)/(i+da)*sum(c(k,obst).*k'.*a(k-1,i+da-1,obst) ) -sum(c(k,obst).*a(k,i+da,obst));
%         end
%         if ~isempty(r)
%             extra = extra -ft(obst,r,1)*a(r-1,i+da,obst);
%         end
%         a(1,i,obst) = (a(1,1,obst)/(i+da)*(ft(obst,2,i+da) + a(1,i,obst)  ) ...
%             + extra  + c(i+da,other) -ft(obst,1,i+da+1) -tmp)/denom; %Minus sign is taken into sum
% %         a(1,i,obst) = (a(1,1,obst)/(i+da)*(ft(obst,2,i+da) + a(1,i,obst) + sum(c(k,obst).*k'.*a(k-1,i+da-1,obst) ) ) ...
% %             + extra -sum(c(k,obst).*a(k,i+da,obst)) + c(i+da,other) -ft(obst,r,1)*a(r-1,i+da,obst) -ft(obst,1,i+da+1) -tmp)/denom;
%     end
%     for obst = 1:J
%         for ii = 2:i
%             for l = ii:i %maxOrder
%                 k = (ii-1):(l-1);
%                 a(ii,l,obst) = sum(a(ii-1,k,obst).*a(1,l-k,obst) );
%             end
%         end
%         c(i+da, obst) = 0;
%         for l = 1:i-1
%             n = 3:(i-l+2);
%             c(i+da, obst) = c(i+da, obst) + (n-1).*ft(obst,n,l)*a(n-2,i-l,obst);
%         end
%         k = 2:i+da-1;
%         c(i+da, obst) = -1/i/a(i+da-1,i+da-1,obst)*(ft(obst,2,i+da) + c(i+da, obst) + k*(c(k,obst).*a(k-1,i+da-1,obst)) );
%     end
%     
% %     na = norm(a(1:i,:) );
%     nc = norm(c(1:i+da,:) );
%     na = norm(squeeze(a(1,1:i,:)) );
%     if (nc == inf) || isnan(nc) || (nc == 0) || (na == inf) || isnan(na) || (na == 0)
%         error(['Incorrect results at i = ' num2str(i)]);
%     end
% end
    
end
