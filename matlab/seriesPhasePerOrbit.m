% Compute the series expansion of the phase of the periodic orbit mode
% Input
%	par       - The par structure from getObst, containing an 'obsts' field in which the J obstacles are ordered according to the periodic orbit. When J > 2,
%       this function should be called with all collections of J-1 obstacles, J-2 and so on till J-J+2 to obtain the phases of all
%       periodic orbits in the scene. The order of the obstacles is important as well.
%   maxOrder  - The maximal order of the series expansion
% Output
%   taus  - The J parameters on the respective obstacles giving the periodic orbit
%   c     - The expansion coefficients of the phase
%   a     - The expansion coefficients of the stationary point
%   ft    - Series expansion coefficients of the distance
function [taus, c, a, ft] = seriesPhasePerOrbit(par, maxOrder)

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
% taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'newton', 'SpecifyObjectiveGradient', true, ...
%     'HessianFcn','objective', 'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000));
taus = fminunc(@distWHess, taus, optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, ...
    'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000, 'Display', 'off'));
% [dt, g, H] = distWHess(taus) % check optimality

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
        gob = 2*transpose(difGam)*gam(:,2)/curDis;
        g(ob) = g(ob) + gob;
        gprev = -2*transpose(difGam)*prGam(:,2)/curDis;
        g(prev) = g(prev) + gprev;
        
        H(ob,ob) = H(ob,ob) +(sum(gam(:,3).^2) +transpose(difGam)*gam(:,2).^2 -g(ob)/2/curDis)/curDis;
        H(ob,prev) = H(ob,prev) +(sum(gam(:,2).*prGam(:,2)) - gprev*gob*0.9 )/curDis;
        H(prev,ob) = H(prev,ob) +(sum(gam(:,2).*prGam(:,2)) - gob*gprev )/curDis;
        H(prev,prev) = H(prev,prev)  -(sum(prGam(:,3).^2) +transpose(difGam)*prGam(:,2).^2*1.1 + g(prev)/2/curDis)/curDis;
        prev = ob;
        prGam = gam;
    end
end

binom = @(x,n) prod(x:-1:(x-n+1) )/prod(1:n)*0^((n<0) + abs(rem(n,1)) ); % Not to be used vectorized (and n should be a positive integer)

% Actually, we only need (1, maxOrder, maxOrder) if J = 2 but this would slightly complicate the implementation.
ft = nan(J, maxOrder+1, maxOrder+1);
for obst = 1:J
    % Take next obstacle as in makeV1.m
    other = obst +1 -(obst == length(par.obsts))*length(par.obsts);
    Gama = par.obsts(obst).serpar(taus(obst), maxOrder+1);
    Gamb = par.obsts(other).serpar(taus(other), maxOrder+1);
    Lambda = nan(maxOrder+1, maxOrder+1, 2);
    Lambda(1,1,:) = Gamb(:,1).^2 - 2*Gama(:,1).*Gamb(:,1) + Gama(:,1).^2;
    ft(obst,1,1) = sqrt(sum(Lambda(1,1,:)));
    for l = 2:maxOrder+1
        Lambda(1,l,:) = sum(Gamb(:,1:l).*Gamb(:, l:-1:1), 2) -2*Gama(:,1).*Gamb(:,l);
    end
    for j = 2:maxOrder+1
        Lambda(j,1,:) = sum(Gama(:,1:j).*Gama(:, j:-1:1), 2) -2*Gama(:,j).*Gamb(:,1);
    end
    for l = 2:maxOrder+1
        for j = 2:maxOrder+1
            Lambda(j,l,:) = -2*Gama(:,j).*Gamb(:,l);
        end
    end
    zt = nan(2*maxOrder, maxOrder+1, maxOrder+1);
    % Avoid wrong conversion of dimensions by this loop
    for i = 1:maxOrder+1
        zt(1,i,:) = (Lambda(i,:,1) + Lambda(i,:,2))/ft(obst,1,1)^2;
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
    for k = 1:maxOrder+1
        for i = (1+(k==1)):maxOrder+1
            ft(obst, i,k) = zt(1,i,k)/2;
            for m = 2:k-2+i
                ft(obst,i,k) = ft(obst,i,k) + binom(1/2,m)*zt(m,i,k);
            end
            ft(obst,i,k) = ft(obst,1,1)*ft(obst,i,k);
        end
    end
end

a = nan(J, maxOrder, maxOrder); % a(:,1,maxOrder) should not be used when J==2
c = nan(J, maxOrder);

if J == 2
    % No tolerance 20*eps but 1e-10 due to loss of digits in objective fct in fmincon for quasi-Newton
    % with exact parametrisation without FFT. Even 1e-6*d(1) for trust-region which uses the analytic Hessian.
    tol = ft(1,1,1)*1e-7;
    if (abs(ft(1, 2, 1)) > tol) || (abs(ft(1, 1, 2)) > tol) 
        error('Distance should be minimal.');
    elseif (abs(ft(1, 3, 1)) < tol) && (abs(ft(1, 1, 3)) < tol) && (abs(ft(1, 2,2)) < tol)
        warning('Cubic equation is not implemented.');
    elseif (abs(ft(1, 3, 1)) < tol) || (abs(ft(1, 1, 3)) < tol) || (abs(ft(1, 2,2)) < tol)
        error('No limit cycle exists.');
    elseif (ft(1,3,1) < 0) || (ft(1,1,3) < 0)
        error('Distance should be minimal.')
    end
    a(1,1,1) = -2*ft(1,1,3)/ft(1,2,2) + sign(ft(1,2,2))*sqrt((2*ft(1,1,3)/ft(1,2,2))^2 - ft(1,1,3)/ft(1,3,1));
    a(2,1,1) = -ft(1,2,2)/(a(1,1,1)*ft(1,2,2) +4*ft(1,1,3));
    c(1,1) = 0;
    c(2,1) = 0;
    c(1,2) = -ft(1,2,2)/2/a(1,1,1) - ft(1,3,1);
    c(2,2) = -ft(1,2,2)/2/a(2,1,1) - ft(1,1,3);
    for i = 3:maxOrder
        for l = (i-1):maxOrder 
            k = (l-i+2):(l-1);
            a(:,l-i+3,l) = sum(a(:,l-i+2,k).*a(:,1,l-k), 3);
        end
        M = zeros(4,4);
        M(1,1) = a(1,i,i);
        M(1,2) = -1; 
        % Both M(3,1) and M(4,2) should be zero by construction
        M(2,2) = a(2,i,i);
        M(2,1) = -1;
        M(3,1) = i*a(1,i-1,i-1);
        M(3,3) = 2*(ft(1,3,1) + c(1,2));
        M(4,2) = i*a(2,i-1,i-1);
        M(4,4) = 2*(ft(1,1,3) + c(2,2));
        rhs = -[(ft(1,1,i+1)+ sum(c(1, 3:i-1).*a(1, 3:i-1, i))); ...
            (ft(1,i+1,1) + sum(c(2, 3:i-1).*a(2, 3:i-1, i))); ...
            (ft(1,2,i) +sum((3:i-1).*c(1, 3:i-1).*a(1, 2:i-2,i-1))) ; ...
            (ft(1,i,2) +sum((3:i-1).*c(2, 3:i-1).*a(2, 2:i-2,i-1)))];
        rhs(1) = rhs(1) - (ft(1,3,1) + c(1,2))*sum(a(1, 1, 2:i-2).*a(1, 1, i-2:-1:2)) - sum(ft(1,4:i+1,1).*a(1, 3:i,i));
        rhs(2) = rhs(2) - (ft(1,1,3) + c(2,2))*sum(a(2, 1, 2:i-2).*a(2, 1, i-2:-1:2)) - sum(squeeze(ft(1,1,4:i+1)).*a(2, 3:i, i)');
        for l=1:i-1
            rhs(1) = rhs(1) - sum(a(1, l, l:(i-1-(l==1))).*ft(1, l+1, (i+1-l):-1:(2+(l==1)) ) );
            rhs(2) = rhs(2) - ft(1,(i+1-l):-1:(2+(l==1)), l+1)*squeeze(a(2, l, l:(i-1-(l==1))));
            rhs(3) = rhs(3) -sum(ft(1, (3+(l==1)):(i-l+2), l).*((2+(l==1)):(i-l+1)).*a(1, (1+(l==1)):(i-l), i-l));
            rhs(4) = rhs(4) -(((2+(l==1)):(i-l+1)).*a(2, (1+(l==1)):(i-l), i-l))*squeeze(ft(1,l,(3+(l==1)):(i-l+2)));
        end
        s = M\rhs;
        c(:,i) = s(1:2);
        a(:,1,i-1) = s(3:4);
    end
    % Could return now but continue to check whether this corresponds to the general case J >= 2
end

for obst = 1:J
    c(obst,1) = -ft(obst,2,1);
end
if (J == 2) && (norm(c(:,1) ) > tol)
    error('taus does not give the minimum of the distance.')
elseif (J ~= 2) && (norm(ft(:,1,2) -[c(2:end,1); c(1,1)]) > norm(c(:,1))*eps^(3/7))
    error('taus does not give the minimum of the distance.')
end

    % op = [omega_{1,2}, omega_{2,2}, ... ; psi_{1,1}, psi_{2,1}, ...]
    % ca = [c_{1,2}, c_{2,2}, ... ; a_{1,1,1}, a_{2,1,1}... ]
    % Actually, \omega_{j,2} = 0 = [c_{j,2} + f_{j,3,1} ] (a_{j,1,1})^2 + f_{j,2,2} a_{j,1,1}  - c_{j+1,2}  + f_{j,1,3}, 
    % but we get [c_{j,2} + f_{j,3,1} ] = -f_{j,2,2}/2/a_{j,1,1} from psi_{j,1} to simplify the formulae to 
    % \omega_{j,2} = 0 = f_{j,2,2} a_{j,1,1}/2  - c_{j+1,2}  + f_{j,1,3}, 
    function op = omPsi(ca)
        op = [[-ca(1,2:end), -ca(1,1)]; 2*ca(1,:).*ca(2,:)];
        for ob = 1:size(ca,2)
            op(1,ob) = op(1,ob) +ft(ob,1,3) + ft(ob,2,2)/2*ca(2,ob);
            op(2,ob) = op(2,ob)  + ft(ob,2,2) +2*ft(ob,3,1)*ca(2,ob);
        end
    end

ca = fsolve(@omPsi, [+transpose(ft(:,3,1)); zeros(1,J)], optimoptions('fsolve', ...
    'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000, 'Display', 'off') );
if (norm(c(:,2)' -ca(1,:)) > eps^(1/2)*norm(c(:,2)) ) || (norm(a(:,1,1) -ca(2,:)') > eps^(1/2)*norm(ca(2,:)) )
    error('Different result'); % if c not set then nan so tests fail
end
c(:,2) = ca(1,:);
a(:,1,1) = ca(2,:);
if (J ~= 2) && any(a(:,1,1) > abs(ft(:,1,2)./ft(:,2,1)) )
    warning('\chi does not contract around the periodic orbit');
elseif any(c(:,2) < eps^(1/2))
    warning('c_{j,2} is zero or negative so \tau_{j+1}^* is not close to a minimum');
end

for i = 3:maxOrder
    for l = (i-1):maxOrder
        k = (l-i+2):(l-1);
        a(:,l-i+3,l) = sum(a(:,l-i+2,k).*a(:,1,l-k), 3);
    end
    M = zeros(2*J, 2*J);
    rhs = zeros(2*J,1);
    for obst = 1:J
        other = obst +1 -(obst == J)*J;
        M(obst, obst) = a(obst, i,i);
        M(obst, other) = -1;
        % a_{j,1,i-1} should not influence omega_{j,i} because of omega_{j,2}
        M(J+obst, obst) = i*a(obst, i-1, i-1);
        M(J+obst, J+obst) = 2*(ft(obst,3,1) + c(obst,2));
        
        rhs(obst) = -sum(c(obst,3:i-1).*a(obst,3:i-1,i)) - sum(ft(obst,4:i+1,1).*a(obst, 3:i,i)) ...
           -ft(obst,1,i+1) - (ft(obst,3,1) + c(obst,2))*sum(a(obst, 1, 2:i-2).*a(obst, 1, i-2:-1:2));
        rhs(J+obst) = -ft(obst,2,i) - sum((3:i-1).*c(obst, 3:i-1).*a(obst, 2:i-2, i-1));
        
        for l=1:i-1
            rhs(obst) = rhs(obst) - sum(a(obst, l, l:(i-1-(l==1))).*ft(obst, l+1, (i+1-l):-1:(2+(l==1)) ) );
            rhs(J+obst) = rhs(J+obst) -sum(ft(obst, (3+(l==1)):(i-l+2),l).*((2+(l==1)):(i-l+1)).*a(obst, (1+(l==1)):(i-l), i-l));
        end
    end
    s = M\rhs; 
    if (J == 2) && ( (norm(c(:,i) -s(1:2)) > eps^(1/2)*norm(c(:,i)) ) || (norm(a(:,1,i-1) -s(3:4)) > eps^(1/2)*norm(s(3:4))) )
        error('Different result');
    end
    c(:,i) = s(1:J);
    a(:,1,i-1) = s(J+1:end);
end

    
end
