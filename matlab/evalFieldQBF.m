% Evaluate the scattered field.
% Input
%   par		  - The structure containing k,N,par,gradnorm,corners,t,colltau
%   x		  - Point where to evaluate
%   sol		  - Solution vector
%   comprSol 	  - 1 if the solution is compressed (for possibly adding the phase), else 0
% Output
%   z		  - The scattered field
function z = evalFieldQBF(par,x,sol,comprSol)

if isfield(par,'obsts')
    par.dbf = par.obsts(1).dbf; % Assume the same degree of basis functions to avoid refactoring.
end

if par.dbf == 1
    if isfield(par, 'quadR') && isfield(par, 'fco')
        % Use a simple Riemann sum with the given number of points per collocation points interval to compute the integral.
        % Using the middle sum avoids points of nonanalyticity.
        if comprSol
            qbf_x = linspace(-1, 1, par.quadR*2+1);
        else
            qbf_x = linspace(-1, 1, round(par.quadR*length(par.fco)/par.N)*2+1);
        end
        step = qbf_x(2)-qbf_x(1);
        qbf_x = qbf_x(1:end-1) + step/2;
        qbf_w = (1-abs(qbf_x))/sum(1-abs(qbf_x));
        istep = round(1/step);
    elseif isfield(par, 'quadR') && isfield(par, 'obsts') && isfield(par.obsts(1), 'fco')
        if comprSol
            qbf_x = linspace(-1, 1, par.quadR*2+1);
        else
            qbf_x = linspace(-1, 1, round(par.quadR*length(par.obsts(1).fco)/par.obsts(1).N)*2+1);
        end
        step = qbf_x(2)-qbf_x(1);
        qbf_x = qbf_x(1:end-1) + step/2;
        qbf_w = (1-abs(qbf_x))/sum(1-abs(qbf_x));
        istep = round(1/step);
    else
        qbf_w = [   0.01666666666667   0.26666666666667   0.43333333333333   0.26666666666667   0.01666666666667];
        step = 1/2; istep = 1/step;
    end
    dbf = 1;
elseif par.dbf == 3  % Weights computed in iepack/basis/make_bfinfo_periodic_spline_equi.m using quad_sf with lo_d = sqrt(2)/16*[1 4 6 4 1]
	qbf_w = [  0.000022045855379   0.009876543209876   0.084744268077601  0.238447971781305   0.333818342151675  ...
		0.238447971781305  0.084744268077601   0.009876543209877   0.000022045855379];
	step = 4/(length(qbf_w)-1);
	istep = round(1/step);
    dbf = 3;
end
if isfield(par,'obsts') && isfield(par.obsts(1), 'phase') && comprSol
	z = 0;
	for obst = 1:length(par.obsts)
		tau = ((0:par.obsts(obst).fN*istep-1)*step -dbf)/par.obsts(obst).fN;
		% Add points to the right using periodicity.
		tau2 = [tau tau(1:dbf*istep + 1)];
        
		kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.obsts(obst).par(tau2) ).^2, 1)));
		gnvals = par.obsts(obst).gradnorm(tau2).*exp(1i*par.k*par.obsts(obst).phase(tau2) );
		zt = 0;
		for i=1:par.obsts(obst).fN
			zt = zt + sol(i+par.fr(1,obst)-1) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
				.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
		end
		z = z + zt/par.obsts(obst).fN;
	end
elseif isfield(par,'obsts') % Ordinary multiple scattering
	z = 0;
    if length(sol) ~= par.N, error('sizes should be compatible for ordinary multiple scattering'); end
	for obst = 1:length(par.obsts)
		tau = ((0:par.obsts(obst).N*istep-1)*step -dbf)/par.obsts(obst).N;
		% Add points to the right using periodicity.
		tau2 = [tau tau(1:dbf*istep + 1)];
        
		kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.obsts(obst).par(tau2) ).^2, 1)));
		gnvals = par.obsts(obst).gradnorm(tau2);
		zt = 0;
		for i=1:par.obsts(obst).N
			zt = zt + sol(i+par.r(1,obst)-1) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
				.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
		end
		z = z + zt/par.obsts(obst).N;
	end
elseif comprSol && isfield(par,'difase')
    fN = length(par.fco);
	tau = ((0:fN*istep-1)*step -dbf)/fN;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
	kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.par(tau2) ).^2, 1)));
	gnvals = par.gradnorm(tau2).*exp(2i*pi*par.k*interp1(par.Tase,par.difase,tau2, 'spline', 'extrap'));
	z = 0;
	for i=1:fN
		z = z + sol(i) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
			.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
	end
	z = z/fN;
elseif comprSol && isfield(par,'phase')
    fN = length(par.fco);
	tau = ((0:fN*istep-1)*step -dbf)/fN;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
	kernelvals = 1i/4*besselh(0,1,par.k*sqrt(sum((repmat(x,1,size(tau2,2))-par.par(tau2) ).^2, 1)));
	gnvals = par.gradnorm(tau2).*exp(1i*par.k*par.phase(tau2));
	z = 0;
	for i=1:fN
		z = z + sol(i) * sum(qbf_w .* kernelvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) ...
			.* gnvals((i-1)*istep+1:(i-1)*istep+length(qbf_w)));
	end
	z = z/fN;
else
	tau = ((0:par.N*istep-1)*step -dbf)/par.N;
	% Add points to the right using periodicity.
	tau2 = [tau tau(1:dbf*istep + 1)];
	
    kgVals = kernelGrad(par, tau2, x);
	z = 0;
	for i=1:par.N
		z = z + sol(i) * sum(qbf_w .* kgVals((i-1)*istep+1:(i-1)*istep+length(qbf_w)) );
	end
	z = z/par.N;
end

end

% Evaluate the kernel for the single or double layer potential times the gradient of the parametrization.
function ker = kernelGrad(par, tau, x)

x = repmat(x,1,length(tau));
if isfield(par, 'secondKind') % dlp
    xy = x-par.par(tau);
    res = sqrt(sum(xy.^2, 1));
    z = 0*res;
    
    I1 = find(abs(res) > 1e-12); % I1 can be empty
    z(I1) = sum(par.normal(tau(I1)) .*xy(:,I1)) .*1i/4 .*par.k .*besselh(1, 1, par.k*res(I1)) ./res(I1);
    
    I2 = find(abs(res) <= 1e-12); % I2 can be empty
    g = par.grad(tau(I2));
    dg = par.derGrad(tau(I2));
    z(I2) = (dg(1,:).*g(2,:) -g(1,:).*dg(2,:))./(4*(g(1,:).^2 +g(2,:).^2).^(3/2) *pi);
    ker = z.*par.gradnorm(tau);
else % slp
    ker = 1i/4.*besselh(0, 1, par.k*sqrt(sum((x-par.par(tau) ).^2, 1)) ).*par.gradnorm(tau);
end
end