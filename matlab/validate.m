% Validate the BEM-matrices. 
% Input
%   A1		- BEM matrix without compression
%   A2		- BEM matrix with compression
%   par		- The structure containing k,N,par,gradnorm,corners,t,colltau, bc, dbf, ppw and xi
%   v		- Validation data: 2 refers to A and A2=\tilde{A}, 2+mti to A\b, A2\b and the solutions of the iterative solver.
%		       conds(:,1:2) contains the condition numbers, errBCavm(:,1:2+mti) contains 1-norm errors on BC on random places given by v.taus, 
%		       errTrueF(:,1:2+mti) contains errors with respect to the true solution of the unit circle on v.rtests and v.taus, 
%		       nnz & perc(:,1) have the (fraction of) nonzeros, errSol(:,1:2+mti) the error wrt c1, errBCcol(:,1:2+mti) the error on the BC on 
%		       the coll pts (= A1*c-b), compresErr(:,1:2+mti) the error made by compressing A*c1-b, timeSol(:,1:2+mti) the time to compute x, 
%		       nbIter(:,1:nbGm) has the number of iterations of the iterative solver, timeA(:,1:2) the times taken to compute A previously and 
%		       field(ny,nx) the total field value at [v.xs(nx),v.ys(ny)]
%   ix		- The row-index where to update the validation data
% Output
%   v		- Updated validation data
function v = validate(A1,A2,par,v,ix)
if isfield(v,'nnz')
	v.nnz(ix,1) = nnz(A1);
	v.nnz(ix,2) = nnz(A2);
end
if isfield(v,'perc')
	v.perc(ix,1) = nnz(A1)/numel(A1);
	v.perc(ix,2) = nnz(A2)/numel(A2);
end
if isfield(v,'conds')
    if (size(A1,1) > 3e3)
        display([datestr(now) ': Starting computing condition numbers ...']);
    end
    % We can also use condest for large matrices to reduce time but the results are sometimes an order of magnitude different.
	v.conds(ix,1) = cond(A1);
	v.conds(ix,2) = cond(A2);
    if (size(A1,1) > 3e3)
        display(['... Ended computing condition numbers on ' datestr(now)]);
    end
end
if isfield(par,'obsts') % Multiple scattering objects
	b = zeros(par.N,1);
	for obst = 1:length(par.obsts)
		b(par.r(1,obst):par.r(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
	end
else
	b = par.bc(par.k,par.par(par.colltau));
end
start = toc; % Avoid restarting the clock with tic because the calling function might already be using it.
c1 = A1\b;
if isfield(v,'timeSol'), v.timeSol(ix,1) = toc-start; end
if isfield(v,'errBCcol'), v.errBCcol(ix,1) = norm(A1*c1-b)/norm(b); end % Not that meaningful but for completeness
if isfield(v,'compresErr'), v.compresErr(ix,1) = norm(A1*c1-b)/norm(b); end % Not that meaningful but for completeness

exSolField = NaN;
if isfield(v,'errTrueF')
    % Use a slight adjustment of the functions in withoutDiscret.m
    n = 0;
    a = par.par(0);
    a = a(1);
    osv = ones(size(v.rtests));
    prevTerm = osv;
    term = -besselj(n,par.k*a*ones(size(v.rtests)))./besselh(n,1,par.k*a*ones(size(v.rtests))).*besselh(n,1,par.k*v.rtests);
    exSolField = term;
    n = n+1;
    while (norm(term) +norm(prevTerm) > 1e-5) && (n < 1.7*par.k+50)
        prevTerm = term;
        term = 2*(1i)^n*(-besselj(n,par.k*a*osv)./besselh(n,1,par.k*a*osv).*besselh(n,1,par.k*v.rtests) ).*cos(2*pi*n*v.taus);
        exSolField = exSolField+term;
        n = n+1;
    end
    if (norm(term) +norm(prevTerm) > 1e-5) || isnan(norm(exSolField) )
        warning([num2str(n) ' = n, exact solution in the field did not converge at k = ' num2str(par.k) ]);
    end
    n = 0;
    term = besselj(n,par.k*a*ones(par.N,1))./besselh(n,1,par.k*a*ones(par.N,1)).*par.k.*...
		(besselh(n-1,1,par.k*a*ones(par.N,1)) -n/(par.k*a)*besselh(n,1,par.k*a*ones(par.N,1) ));
    % c3 is an interpolation of the solution.
    c3 = term -1i*par.k*cos(2*pi*par.colltau').*exp(1i*par.k*a*cos(2*pi*par.colltau'));
    prevTerm = ones(par.N,1);
    n = n+1;
    while (norm(term) + norm(prevTerm) > 1e-5) && (n < 1.7*par.k+50)
        prevTerm = term;
        term = 2*(1i)^n*besselj(n,par.k*a*ones(par.N,1))./besselh(n,1,par.k*a*ones(par.N,1)).*par.k.*...
		(besselh(n-1,1,par.k*a*ones(par.N,1)) -n/(par.k*a)*besselh(n,1,par.k*a*ones(par.N,1) )).*cos(2*pi*n*par.colltau');
        c3 = c3 +term;
        n = n+1;
    end % c3 should be close to c1, but A1 is probably ill-conditioned and we use c3 only here.
    if (norm(term) + norm(prevTerm) > 1e-4) || isnan(norm(c3) )
        warning([num2str(n) ' = n, c3 did not converge at k = ' num2str(par.k) ]);
    end
    for av=1:v.avm
        point = [v.rtests(av)*cos(2*pi*v.taus(av)); v.rtests(av)*sin(2*pi*v.taus(av))];
        v.errTrueF(ix,1) = v.errTrueF(ix,1) + abs(exSolField(av) - evalFieldQBF(par, point, c1, 0) )/v.avm;
        v.errTrueF(ix,2) = v.errTrueF(ix,2) + abs(exSolField(av) - evalFieldQBF(par, point, c3, 0) )/v.avm;
    end
end

start = toc;
if isfield(par, 'fco') % One can add this field to par when the phase is known and N can be reduced
    b2 = par.bc(par.k,par.par(par.fco));
    if 1
        c2 = A2\b2;
    else
       c2 = -par.k*2i.*cos(2*pi*par.fco').*(par.fco' > 0.25).*(par.fco' < 0.75); % plane wave on a circle
    end
    figure; plot(par.fco, [real(c2) imag(c2)]); legend('Re(c2) using fco','Im(c2) using fco')
elseif isfield(par, 'obsts') && isfield(par.obsts(1), 'fco')
    b2t = nan*zeros(size(A2,1),1);
    for obst = 1:length(par.obsts)
%         b2t(2*(par.r(1,obst):par.r(2,obst))-1) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
%         b2t(2*(par.r(1,obst):par.r(2,obst))) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau + par.obsts(obst).hs));
        % Above for two phases, below for one
%         b2t(par.r(1,obst):par.f(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).colltau));
        b2t(par.fr(1,obst):par.fr(2,obst)) = par.bc(par.k,par.obsts(obst).par(par.obsts(obst).fco));
    end
    c2 = A2\b2t;
    warning('b2t might not be the expected rhs.');
else
    b2 = b;
    c2 = A2\b;
end
if isfield(v, 'timeSol'), v.timeSol(ix,2) = toc-start; end
if isfield(v, 'errSol')
    if isfield(par, 'phase')
%         v.errSol(ix,2) = norm(interp1(par.fco, c2, par.colltau').*exp(1i*par.k*par.phase(par.colltau) )-c1)/norm(c1);  
%         v.errSol(ix,2) = norm(interp1(par.fco, c2, par.colltau').*transpose(exp(1i*par.k*par.phase(par.colltau) ))-c1)/norm(c1);
        v.errSol(ix,2) = norm(interp1([par.fco (par.fco(1)+1)], [c2; c2(1)], par.colltau').*transpose(exp(1i*par.k*par.phase(par.colltau) ))-c1)/norm(c1);
    elseif isfield(par, 'obsts') && isfield(par.obsts(1), 'phase')
        c2ip = nan*c1;
        for obst = 1:length(par.obsts)
            c2ip(par.r(1,obst):par.r(2,obst)) = interp1([par.obsts(obst).fco (par.obsts(obst).fco(1)+1)], [c2(par.fr(1,obst):par.fr(2,obst) ); ...
                c2(par.fr(1,obst) )], par.obsts(obst).colltau').*transpose(exp(1i*par.k*par.obsts(obst).phase(par.obsts(obst).colltau) ));
            debug = 1; % any(isnan(c2ip(par.r(1,obst):par.r(2,obst)) ))
        end
        v.errSol(ix,2) = norm(c2ip-c1)/norm(c1);
    else
        v.errSol(ix,2) = norm(c2-c1)/norm(c1); 
    end
end
if isfield(v, 'errBCcol'), v.errBCcol(ix,2) = norm(A1*c2-b)/norm(b); end
if isfield(v, 'compresErr'), v.compresErr(ix,2) = norm(A2*c1-b)/norm(b); end

if isfield(v, 'errBCavm')
    v.errBCavm(ix,:) = 0;
    for av = 1:v.avm
        if isfield(par,'obsts'), x = par.obsts(mod(round(v.taus(av)/eps),length(par.obsts))+1).par(v.taus(av));
        else x=par.par(v.taus(av));
        end
        v.errBCavm(ix,1) = v.errBCavm(ix,1) + abs(par.bc(par.k,x) - evalFieldQBF(par,x,c1,0) )/v.avm;
        v.errBCavm(ix,2) = v.errBCavm(ix,2) + abs(par.bc(par.k,x) - evalFieldQBF(par,x,c2,1) )/v.avm;
    end
end

if isfield(v, 'errInt')
    v.intXs = NaN*zeros(2,v.avm); % We could reuse these evaluation points for each k but would have to reset for new objects.
    mnr = 200;
    ts = linspace(0,1,mnr);
    % Take a line between two random points on the boundary and ensure the connecting line is not outside a nonconvex obstacle.
    warning 'off'; % Avoid warnings of ill-conditioned system when computing intersections.
    for av = 1:v.avm
        if isfield(par,'obsts') % For multiple scattering obstacles, take a random obstacle given by the (fixed) v.taus
            cpar = par.obsts(mod(round(v.taus(av)/eps),length(par.obsts))+1);
        else cpar=par;
        end
        parametr = [cpar.par(ts), cpar.par(0)];
        mi = min(parametr,[],2);
        ma = max(parametr,[],2);
        outbb = ma + 0.1*(ma-mi);
        out = 1;
        while out
            v.intXs(:,av) = mi + diag(rand(2,1))*(ma-mi);
            % Make a line to the exterior and if the number of jumps over the boundary is odd then we are inside.
            for mni = 1:mnr
                sol = [ (parametr(:,mni+1)-parametr(:,mni)), (v.intXs(:,av)-outbb)]\(v.intXs(:,av)-parametr(:,mni));
                if ( (min(sol) > 0) && (max(sol) < 1) )
                    out = ~out;
                end
            end
        end
    end
    warning 'on';
    v.errInt(ix,:) = 0;
    for av = 1:v.avm
        x = v.intXs(:,av);
        v.errInt(ix,1) = v.errInt(ix,1) + abs(par.bc(par.k,x) - evalFieldQBF(par,x,c1,0) )/v.avm;
        v.errInt(ix,2) = v.errInt(ix,2) + abs(par.bc(par.k,x) - evalFieldQBF(par,x,c2,1) )/v.avm;
    end
end

tol = 1e-5; % Or take a (GMRES) tolerance v.errSol(ix,2)/10 or decreasing as 1e-2/par.k^2 or 1/cond(A1).
itproc = @gmres; % We can choose gmres, bicgstabl, bicgstab, ....
if isfield(v,'mti'), for ti = 1:v.mti
    start = toc;
    % A1\b, A2\b, A1\b with precond A2, A1\b with x0=A2\b
    if ti == 1, [c,flag,relres,iter,resvec] = itproc(A1,b,par.N,tol,par.N);
    elseif ti == 2, [c,flag,relres,iter,resvec] = itproc(A2,b2,par.N,tol,par.N);
    elseif ti == 3, [c,flag,relres,iter,resvec] = itproc(A1,b,par.N,tol,par.N,A2);
    elseif ti == 4, [c,flag,relres,iter,resvec] = itproc(A1,b,par.N,tol,par.N,[],[],c2);
	else error('Too many types asked for the iterative solver'); end;
    if isfield(v,'timeSol'), v.timeSol(ix,2+ti) = toc-start; end
    if (flag ~= 0) || (relres > tol)
        warning([num2str(ti) ' = ti, GMRES wrong at ix = ' num2str(ix) ' with flag = ' num2str(flag) ...
            ', relres = ' num2str(relres) ' and tol = ' num2str(tol)]);
    end
    if isfield(v,'nbIter'), v.nbIter(ix,ti) = iter(2); end
    if isfield(v,'errSol'); v.errSol(ix,2+ti) = norm(c-c1)/norm(c1); end
    if isfield(v,'errBCcol'), v.errBCcol(ix,2+ti) = norm(A1*c-b)/norm(b); end
    if isfield(v, 'errBCavm'), for av = 1:v.avm, 
		if isfield(par,'obsts'), x = par.obsts(mod(round(v.taus(av)/eps),length(par.obsts))+1).par(v.taus(av));
		else x=par.par(v.taus(av));
		end
		v.errBCavm(ix,2+ti) = v.errBCavm(ix,2+ti) + abs(par.bc(par.k,x) - evalFieldQBF(par,x,c,ti==2) )/v.avm;
    end, end
    if isfield(v, 'errInt'), for av = 1:v.avm, 
		v.errInt(ix,2+ti) = v.errInt(ix,2+ti) + abs(par.bc(par.k,v.intXs(:,av)) - evalFieldQBF(par,v.intXs(:,av),c,ti==2) )/v.avm;
    end, end
    if isfield(v,'errTrueF') % We cannot have multiple obstacles: this exact solution is for incident x-plane wave on unit circle.
        for av=1:v.avm
            point = [v.rtests(av)*cos(2*pi*v.taus(av)); v.rtests(av)*sin(2*pi*v.taus(av))];
            v.errTrueF(ix,2+ti) = v.errTrueF(ix,2+ti) + abs(exSolField(av) - evalFieldQBF(par,point,c,ti == 2) )/v.avm;
        end
    end
end, end

if isfield(v,'field')
    % Compute the field in a region around the obstacle for plotting.
    nd = size(v.field,1);
    if isfield(par,'obsts')
        v.parametr = cell(length(par.obsts),1);
        mi = [Inf; Inf];
        ma = [-Inf; -Inf];
        for obst = 1:length(par.obsts)
            vp = par.obsts(obst).par(linspace(0,1,round(nd/length(par.obsts))));
            mi = min([mi,vp], [], 2);
            ma = max([ma,vp], [], 2);
            v.parametr{obst} = vp;
        end
    else
        v.parametr = par.par(linspace(0,1,nd));
        mi = min(v.parametr,[],2);
        ma = max(v.parametr,[],2);
    end
    if 1
        siz = 1.1*max(ma-mi);
        v.xs = (mi(1)+ma(1))/2 +linspace(-siz,siz,nd+1);
        v.ys = (mi(2)+ma(2))/2 +linspace(-siz,siz,nd);
    else
        siz = 0.6*(ma-mi);
        v.xs = (mi(1)+ma(1))/2 +linspace(-siz(1), siz(1), nd+1);
        v.ys = (mi(2)+ma(2))/2 +linspace(-siz(2), siz(2), nd);
    end
	for nx = 1:nd+1, for ny = 1:nd
		v.field(ny,nx) = evalFieldQBF(par,[v.xs(nx); v.ys(ny)],c1,0) - par.bc(par.k,[v.xs(nx); v.ys(ny)]); % BC = -Inc field
	end, end
end
