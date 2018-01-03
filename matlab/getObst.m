% Make a parameterisation. Since the wavenumber is variable, it should be added by the calling function as par.k, and then par.N, 
% par.t (knots of the periodic spline) and par.colltau (the collocation points) should be set.
% Input
%	idx	- Index for which obstacle: 0 circle (cubic), 1 circle (linear as are the following), 2 ellipse, 3 near-inclusion, 
%         4 almost convex, 5 two circles, 6 three ellipses, 7 near-inclusion and circle, 8 nonconvex polygon, 9 one selfreflection, 11 two
%         ellipses
%   [rx - the radius in the x-direction for the two ellipses, or the shift for smooth parametrizations based on an FFT]
% Output
%	par - The par structure, containing par (the parametrization), gradnorm, normal, corners, bc (the boundary condition), dbf (degree
%         of the basis functions), ppw (number of points per wavelength), serpar (the series expansion of the parametrization for one t
          % and arbitrary order) and xi (threshold percentage). For multiple scattering obstacles, these are combined in an 'obsts' field.
function par = getObst(idx, rx)
switch idx
	case 0 % Circle with cubic basis functions
		par = struct('par', @(t) repmat([0;0], 1, length(t)) + [0.5*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm', @(t) pi*ones(size(t)), ...
            'normal', @(t) [cos(2*pi*t); sin(2*pi*t)], 'corners', [], 'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'dbf', 3, ...
            'ppw', 6, 'xi', 0.04);
		return
	case 1 % Circle with linear basis functions
		par = struct('par', @(t) repmat([0;0], 1, length(t)) + [0.5*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm', @(t) pi*ones(size(t)), ...
            'normal', @(t) [cos(2*pi*t); sin(2*pi*t)], 'grad', @(t) [-pi*sin(2*pi*t); pi*cos(2*pi*t)] , ...
            'derGrad', @(t) repmat([0;0], 1, length(t)) + [-2*pi^2*cos(2*pi*t); -2*pi^2*sin(2*pi*t)], 'corners', [], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'dbf', 1, 'ppw', 6, 'xi', 0.04);
		return
	case 2 % Ellipse
		par = struct('par', @(t) repmat([0;0], 1, length(t)) + [0.3*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm', ...
			@(t) 2*pi*sqrt(0.3^2+0.5^2)*ones(size(t)), 'corners', [], 'bc', @(k,x) -1*exp(1i*k*(x')*[cos(-pi/4); sin(-pi/4)]), ...
            'dbf', 1, 'ppw', 6, 'xi', 0.04);
		return
	case 3 % Near-inclusion
		x = [0.7454 0.6947 0.5933 0.3836 0.2108 0.1440 0.2546 0.4597 0.5472 0.5219 0.3491 0.1694 0.2131 0.4366 0.6417 0.7200;...
			0.5219 0.6506 0.7646 0.8904 0.8553 0.6857 0.6009 0.6272 0.5395 0.3348 0.3640 0.3319 0.1360 0.1038 0.1974 0.3787];
        dbf = 1; ppw = 10; xi = 0.01;
	case 4 % Near-convex
		x = [ 0.6071  0.6141  0.6302  0.6233  0.5611  0.4965  0.4343  0.3906  0.3836  0.3975  0.4067  0.3952  0.3813  ...
			0.4067  0.4712  0.5426  0.5887  0.6256  0.6141  0.6071 ;   0.5365  0.6272  0.7178  0.8173  0.8874  0.9020  0.8787  ...
			0.8085  0.7120  0.6067  0.5073  0.4079  0.3056  0.2354  0.1827  0.1798  0.2149  0.2909  0.3845  0.4664];
        
        dbf = 1; ppw = 7; xi = 0.04; %xi = 0.025;
	case 5 % Two circles for the article periodic orbit general obstacles with d = 1 = 2r
		par = struct('obsts', [struct('par', @(t) [sin(2*pi*t); cos(2*pi*t)]/2, 'gradnorm',... 
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ... 
            (2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)]/2) ...
			struct('par', @(t) repmat([0; 2], 1, length(t)) + [sin(2*pi*t); -cos(2*pi*t)]/2, 'gradnorm',...
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ...
            (4*((0:mo-1) == 0) -(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2)], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
		return
	case 6 % Three ellipses
		par = struct('obsts', [struct('par', @(t) repmat([0;0], 1, length(t)) + [0.3*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(0.3^2+0.5^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [0.3*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2); ...
            (0 + 0.5*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))] ) ...  %%%%%%%%%%%%%%%%%%%%%%%%
			struct('par', @(t) repmat([-1;1], 1, length(t)) + 0.4*[cos(2*pi*t); sin(2*pi*t)], 'gradnorm',...
			@(t) 0.8*pi*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [ (-1*((0:mo-1) == 0) + 0.4*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ...
            (1*((0:mo-1) == 0) + 0.4*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))]) ...  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			struct('par', @(t) repmat([1;1], 1, length(t)) + [0.4*cos(2*pi*t); 0.2*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(0.4^2+0.2^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [ (1*((0:mo-1) == 0) + 0.4*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ...
            (1*((0:mo-1) == 0) + 0.2*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))] ) ], ...  %%%%%%%%%%%%%%%%%%%%%
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
        return
    case 7 % Near-inclusion and circle
        tmp = getObst(3);
        tmp = rmfield(tmp,'bc');
        tmp = rmfield(tmp,'xi');
        tmp = rmfield(tmp,'grad');
        tmp = rmfield(tmp,'cur');
        tmp = rmfield(tmp,'derGrad');
        tmp.ppw = 10;
        par = struct('obsts', [tmp struct('par', @(t) repmat([-0.5;-0.5], 1, length(t)) + [0.5*cos(2*pi*t); 0.5*sin(2*pi*t)], ...
            'gradnorm', @(t) pi*ones(size(t)), 'normal', @(t) [cos(2*pi*t); sin(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10)],...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 1e-3);
		return
    case 8 % Polygon
        y = [1 -1 0 -1  1; ...
            1 1  0 -1 -1];
        dbf = 1; ppw = 9; xi = 0.04;
    case 9 % One self-reflection
        x = [1 -1 0 -1  1; ...
            1 1  0 -1 -1];
        dbf = 1; ppw = 15; xi = 0.04;
	case 10 % One self-reflection line pieces
		par = struct('par', @(t) onePar, 'gradnorm', @(t) oneGrad(t), 'corners', [], 'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), ...
            'dbf', 1, 'ppw', 15, 'xi', 0.04);
		return
    case 11 % Two ellipses
		par = struct('obsts', [struct('par', @(t) repmat([-0.5;-0.5], 1, length(t)) + [rx*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(rx^2+1/4)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10) ...
			struct('par', @(t) repmat([-0.5;1.5], 1, length(t)) + [rx*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(rx^2+1/4)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10)], ...
        	    'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
		return
	case 12 % Three circles of radius 1/2
		par = struct('obsts', [struct('par', @(t) [sin(2*pi*t); cos(2*pi*t)]/2, 'gradnorm',... 
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ... 
            (2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)]/2) ... %End obst 1
			struct('par', @(t) repmat([0; 4], 1, length(t)) + [sin(2*pi*t); -cos(2*pi*t)]/2, 'gradnorm',...
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ...
            (8*((0:mo-1) == 0) -(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2), ... % End obst 2
            struct('par', @(t) repmat([2*tan(pi/3); 2], 1, length(t)) + [sin(2*pi*t); -cos(2*pi*t)]/2, 'gradnorm',...
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(4*tan(pi/3)*((0:mo-1) == 0)+(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)); ...
            (4*((0:mo-1) == 0)-(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2)], 'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), ...
            'xi', 3e-3);
		return
    case 13 % Ellipse and near-convex
        tmp = getObst(4);
        tmp = rmfield(rmfield(tmp, 'xi'), 'bc');
        tmp = rmfield(rmfield(tmp, 'grad'), 'cur');
        tmp = rmfield(rmfield(tmp, 'normal'), 'derGrad');
        tmp.ppw = 10;
		par = struct('obsts', [tmp  struct('par', @(t) repmat([-0.4;0], 1, length(t)) + [0.15*cos(2*pi*t); 0.25*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(0.15^2+0.25^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [(-0.4*((0:mo-1) == 0) + 0.15*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ... 
            0.25*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)]) ], ...
        	    'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
		return
    case 14 % Near-convex, ellipse and near-inclusion
        tmp = getObst(4);
        tmp = rmfield(rmfield(tmp, 'xi'), 'bc');
        tmp = rmfield(rmfield(tmp, 'grad'), 'cur');
        tmp = rmfield(rmfield(tmp, 'normal'), 'derGrad');
        tmp.ppw = 10;
        tm = getObst(3, [+1; -0.6]);
        tm = rmfield(rmfield(tm, 'xi'), 'bc');
        tm = rmfield(rmfield(tm, 'grad'), 'cur');
        tm = rmfield(rmfield(tm, 'normal'), 'derGrad');
        tm.ppw = 10;
        par = struct('obsts', [tmp  struct('par', @(t) repmat([-0.4;0], 1, length(t)) + [0.15*cos(2*pi*t); 0.25*sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*sqrt(0.15^2+0.25^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [(-0.4*((0:mo-1) == 0) + 0.15*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ...
            0.25*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)])  tm ], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
        return
	case 15 % Convex obstacle with fft
		x = [0.4924    0.2348    0.1600    0.1711    0.3580    0.5810    0.6558    0.6904    0.5935; ...
            0.0934    0.1632    0.4104    0.8142    0.9896    0.9028    0.7104    0.3481    0.1538];
        dbf = 1; ppw = 10; xi = 0.01;
    case 16 % Circle with fft
        x = [cos(2*pi*(0:5)/6); sin(2*pi*(0:5)/6)];
        dbf = 1; ppw = 5; xi = 0.1;
	otherwise
		error(['Unknown index ' num2str(idx)]);
end
% Smooth and polygonal parametrizations, from iepack/boundary/param/par_smooth and par_polygon
if exist('x','var')
    % x are the interpolation points, given anti-clockwise (determining the direction of the outward normal).
    N = size(x,2);
    if exist('rx', 'var')
        x = x + repmat(rx,1,N);
    end
    fx = fft(x(1,:))/N;
    fy = fft(x(2,:))/N;
    if mod(N,2) == 0
        n = [0:N/2, -N/2+1:-1];
    else
        n = [0:floor(N/2), -floor(N/2):-1];
    end
    fxd = 1i*2*pi*(n.*fx);
    fyd = 1i*2*pi*(n.*fy);
    
    par = struct('par', @(t) parS(t), 'gradnorm', @(t) sqrt(sum(gradS(t).^2, 1)), 'serpar', @(t,mo) serparS(t,mo), ...
        'normal', @(t) [1 0; 0 -1]*flipud(gradS(t))./repmat(sqrt(sum(gradS(t).^2, 1)),2,1), ...
        'grad', @(t) gradS(t), 'derGrad', @(t) derGradS(t), 'cur', @(t) curS(t), 'corners', [] );
elseif exist('y','var')
    ny = size(y,2);
    
    par = struct('par', @(t) parPolygon(t), 'gradnorm', @(t) sqrt(sum(gradPolygon(t).^2, 1)), ...
        'normal', @(t) [1 0; 0 -1]*flipud(gradPolygon(t))./repmat(sqrt(sum(gradPolygon(t).^2, 1)),2,1), ...
        'corners', linspace(0,1,ny+1) );
end
par.bc = @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]);
par.dbf = dbf; par.ppw = ppw; par.xi = xi;

function p = parS(t) % parametrisation of a smooth obstacle through interpolation
	p = zeros(2, length(t));
	L = length(t);
	
	A = exp(1i*2*pi*repmat(n',1,L).*repmat(t,N,1));

	if mod(N,2) == 1
		p(1,:) = sum( repmat(fx.', 1, L) .* A, 1);
		p(2,:) = sum( repmat(fy.', 1, L) .* A, 1);
	else
		M = N/2;
		p(1,:) = sum( repmat(fx(1:M).', 1, L) .* A(1:M,:), 1);
		p(2,:) = sum( repmat(fy(1:M).', 1, L) .* A(1:M,:), 1);

		p(1,:) = p(1,:) + sum( repmat(fx(M+2:N).', 1, L) .* A(M+2:N,:), 1);
		p(2,:) = p(2,:) + sum( repmat(fy(M+2:N).', 1, L) .* A(M+2:N,:), 1);
		% add the middle one, a cosine
		p(1,:) = p(1,:) + fx(M+1)*cos(2*pi*n(M+1)*t);
		p(2,:) = p(2,:) + fy(M+1)*cos(2*pi*n(M+1)*t);
	end
	if sum(abs(imag(t))) == 0
		p = real(p);
	end
end

function p = serparS(t, mo) % parametrisation of a smooth obstacle through interpolation
    p = nan(2, length(t), mo);
    L = length(t);
    for ord = 1:mo
        A = (1i*2*pi*repmat(n',1,L)).^(ord-1)./factorial(ord-1).*exp(1i*2*pi*repmat(n',1,L).*repmat(t,N,1));
        
        if mod(N,2) == 1
            p(1,:, ord) = sum( repmat(fx.', 1, L) .* A, 1);
            p(2,:, ord) = sum( repmat(fy.', 1, L) .* A, 1);
        else
            M = N/2;
            p(1,:, ord) = sum( repmat(fx(1:M).', 1, L) .* A(1:M,:), 1);
            p(2,:, ord) = sum( repmat(fy(1:M).', 1, L) .* A(1:M,:), 1);
            
            p(1,:, ord) = p(1,:, ord) + sum( repmat(fx(M+2:N).', 1, L) .* A(M+2:N,:), 1);
            p(2,:, ord) = p(2,:, ord) + sum( repmat(fy(M+2:N).', 1, L) .* A(M+2:N,:), 1);
            % add the middle one, a cosine for ord = 0
            p(1,:, ord) = p(1,:, ord) + (2*pi*n(M+1))^(ord-1)./factorial(ord-1).*fx(M+1)*cos(2*pi*n(M+1)*t +pi/2*(ord-1));
            p(2,:, ord) = p(2,:, ord) + (2*pi*n(M+1))^(ord-1)./factorial(ord-1).*fy(M+1)*cos(2*pi*n(M+1)*t +pi/2*(ord-1));
        end
    end
    if sum(abs(imag(t))) == 0
        p = real(p);
    end
    if L==1
        p = squeeze(p);
    end
    if mo == 1
        p = squeeze(p);
    end
end

function p = curS(t) % Curvature of a smooth obstacle through interpolation
	p = zeros(4, length(t));
	L = length(t);
	
	A = 1i*2*pi*repmat(n',1,L).*exp(1i*2*pi*repmat(n',1,L).*repmat(t,N,1));
    
    if mod(N,2) == 1
        p(1,:) = sum( repmat(fx.', 1, L) .* A, 1);
        p(2,:) = sum( repmat(fy.', 1, L) .* A, 1);
        p(3,:) = sum( repmat(fx.', 1, L).*2i.*pi.*repmat(n',1,L).*A, 1);
        p(4,:) = sum( repmat(fy.', 1, L).*2i.*pi.*repmat(n',1,L).*A, 1);
    else
        M = N/2;
        p(1,:) = sum( repmat(fx(1:M).', 1, L) .* A(1:M,:), 1);
        p(2,:) = sum( repmat(fy(1:M).', 1, L) .* A(1:M,:), 1);
        p(3,:) = sum( repmat(fx(1:M).', 1, L).*2i.*pi.*repmat(n(1:M)',1,L).*A(1:M,:), 1);
        p(4,:) = sum( repmat(fy(1:M).', 1, L).*2i.*pi.*repmat(n(1:M)',1,L).*A(1:M,:), 1);
        
        p(1,:) = p(1,:) + sum( repmat(fx(M+2:N).', 1, L) .* A(M+2:N,:), 1);
        p(2,:) = p(2,:) + sum( repmat(fy(M+2:N).', 1, L) .* A(M+2:N,:), 1);
        p(3,:) = p(3,:) + sum( repmat(fx(M+2:N).', 1, L).*2i.*pi.*repmat(n(M+2:N)',1,L).*A(M+2:N,:), 1);
        p(4,:) = p(4,:) + sum( repmat(fy(M+2:N).', 1, L).*2i.*pi.*repmat(n(M+2:N)',1,L).*A(M+2:N,:), 1);
        
        p(1,:) = p(1,:) - fx(M+1)*2*pi*n(M+1)*sin(2*pi*n(M+1)*t);
        p(2,:) = p(2,:) - fy(M+1)*2*pi*n(M+1)*sin(2*pi*n(M+1)*t);
        p(3,:) = p(3,:) - fx(M+1)*(2*pi*n(M+1))^2*cos(2*pi*n(M+1)*t);
        p(4,:) = p(4,:) - fy(M+1)*(2*pi*n(M+1))^2*cos(2*pi*n(M+1)*t);
    end
    if sum(abs(imag(t))) == 0
        p = real(p);
    end
    p = (p(1,:).*p(4,:) - p(2,:).*p(3,:) )./(p(1,:).^2+p(2,:).^2).^(3/2);
end

function g = gradS(t)
	g = zeros(2, length(t));
	if mod(N,2) == 1
		for j=1:N
			g(1,:) = g(1,:) + fxd(j)*exp(1i*2*pi*n(j)*t);
			g(2,:) = g(2,:) + fyd(j)*exp(1i*2*pi*n(j)*t);
		end
	else
		M = N/2;
		g(1,:) = fxd(1)*exp(1i*2*pi*n(1)*t);
		g(2,:) = fyd(1)*exp(1i*2*pi*n(1)*t);
		for j=1:M-1
			g(1,:) = g(1,:) + fxd(1+j)*exp(1i*2*pi*n(1+j)*t);
			g(1,:) = g(1,:) + fxd(M+1+j)*exp(1i*2*pi*n(M+1+j)*t);
			g(2,:) = g(2,:) + fyd(1+j)*exp(1i*2*pi*n(1+j)*t);
			g(2,:) = g(2,:) + fyd(M+1+j)*exp(1i*2*pi*n(M+1+j)*t);
		end
% Old version: g(1,:) = g(1,:) - imag(fxd(M+1))*sin(2*pi*n(M+1)*t);
		g(1,:) = g(1,:) - fx(M+1)*2*pi*n(M+1)*sin(2*pi*n(M+1)*t);
		g(2,:) = g(2,:) - fy(M+1)*2*pi*n(M+1)*sin(2*pi*n(M+1)*t);
	end
	if sum(abs(imag(t))) == 0
		g = real(g);
	end
end

function x = derGradS(t)
    fxds = fx;
    fyds = fy;
    v = 2;
    while v > 0
        fxds = 1i*2*pi *(n .*fxds);
        fyds = 1i*2*pi *(n .*fyds);
        if mod(N,2) == 0
            fxds(N/2+1) = imag(fxds(N/2+1));
            fyds(N/2+1) = imag(fyds(N/2+1));
        end
        v = v-1;
    end
    x = zeros(2, length(t));
    if mod(N,2) == 1
        for j=1:N
            x(1,:) = x(1,:) + fxds(j)*exp(1i*2*pi*n(j)*t);
            x(2,:) = x(2,:) + fyds(j)*exp(1i*2*pi*n(j)*t);
        end
    else
        M = N/2;
        x(1,:) = fxds(1)*exp(1i*2*pi*n(1)*t);
        x(2,:) = fyds(1)*exp(1i*2*pi*n(1)*t);
        for j=1:M-1
            x(1,:) = x(1,:) + fxds(1+j)*exp(1i*2*pi*n(1+j)*t);
            x(1,:) = x(1,:) + fxds(M+1+j)*exp(1i*2*pi*n(M+1+j)*t);
            x(2,:) = x(2,:) + fyds(1+j)*exp(1i*2*pi*n(1+j)*t);
            x(2,:) = x(2,:) + fyds(M+1+j)*exp(1i*2*pi*n(M+1+j)*t);
        end
        x(1,:) = x(1,:) -fxds(M+1)*cos(2*pi*n(M+1)*t);
        x(2,:) = x(2,:) -fyds(M+1)*cos(2*pi*n(M+1)*t);
    end
    if sum(abs(imag(t))) == 0
        x = real(x);
    end
end



function p = parPolygon(t)
    p = NaN*zeros(2,length(t));
    for ti = 1:length(t)
        l = floor(mod(t(ti),1)*ny); % line between y(:,l) and y(:,l+1)
        if l >= ny-1
            p(:,ti) = y(:,ny) + (y(:,1)-y(:,ny))*(mod(t(ti),1)*ny-ny+1);
        else
            p(:,ti) = y(:,l+1) + (y(:,l+2)-y(:,l+1))*(mod(t(ti),1)*ny-l);
        end
    end
end
function g = gradPolygon(t)
    g = NaN*zeros(2,length(t));
    for ti = 1:length(t)
        l = floor(mod(t(ti),1)*ny); % line between y(:,l) and y(:,l+1)
        if l >= ny-1
            g(:,ti) = (y(:,1)-y(:,ny))*ny;
        else
            g(:,ti) = (y(:,l+2)-y(:,l+1))*ny;
        end
    end
end

function p = onePar(t)
    p = nan*zeros(2, length(t));
    ix = (t >= 0.2) && (t <= 0.45);
    p(:,ix) = 5*t(ix)*[1; -1];
    ix = (t >= 0.55) && (t <= 0.8);
    p(:,ix) = 5*t(ix)*[-1; 1];
end

end % function getObst
