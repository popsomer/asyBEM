% Compute and save V_{j,1}, the first eigenvector of M which represents a full cycle of reflections between obstacles.

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% ks = 2^9; 
ks = 2^7;

% obstacle = 5; % Two circles
% obstacle = 12; % Three circles
% obstacle = 13; % Near-convex and ellipse
obstacle = 14; % Near-convex, ellipse and near-inclusion
swap = 6;

bc = 1;
% bcsh = pi/2; % shift in boundary condition

printtoc = 10;
kl = length(ks);
nbOb = 1; %length(obsts);

avm = 100; % Number of random taus to average BC over
v = struct('avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2),  ...
    'errInt', zeros(nbOb*kl,2), 'timeA', zeros(nbOb*kl,4), 'ks', ks); %, 'field', zeros(70) );

%% Computations
start = now;

for ki = 1:kl
    if swap == 0
        par = getObst(obstacle); % Reset par
    elseif swap == 1 % Three circles switched
        obstacle = '12swap';
%         obstacle = 12; 
        par = getObst(obstacle);
        tmppar = par.obsts(2);
        par.obsts(2) = par.obsts(3);
        par.obsts(3) = tmppar;
    elseif swap == 2 %Two nonsymmetric ellipses
        obstacle = '2nsEll';
        par = struct('obsts', [struct('par', @(t) repmat([-0.5;-1], 1, length(t)) + [0.7*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*sqrt(0.7^2+1/4)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [(-0.5*((0:mo-1) == 0) + 0.7*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ... 
            (-1*((0:mo-1) == 0) +0.5*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))]) ...
            struct('par', @(t) repmat([0.4;0], 1, length(t)) + [0.45*cos(2*pi*t); 0.25*sin(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*sqrt(0.45^2+0.25^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [(0.4*((0:mo-1) == 0) + 0.45*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ... 
            0.25*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)]) ], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
    elseif swap == 3 % Three unequal nonsymmetric circles
        obstacle = '3nsCirc';
		par = struct('obsts', [struct('par', @(t) [sin(2*pi*t); cos(2*pi*t)]/2, 'gradnorm',... 
			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ... 
            (2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)]/2) ... %End obst 1
			struct('par', @(t) repmat([0; 4], 1, length(t)) + 0.35*[sin(2*pi*t); -cos(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*0.35*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [0.35*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ...
            (4*((0:mo-1) == 0) -0.35*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]), ... % End obst 2
            struct('par', @(t) repmat([7.5; 3], 1, length(t)) + 1.1*[sin(2*pi*t); -cos(2*pi*t)], 'gradnorm',...
			@(t) 2*pi*1.1*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
            'serpar', @(t,mo) [(7.5*((0:mo-1) == 0) +1.1*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)); ...
            (3*((0:mo-1) == 0) -1.1*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2)], 'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), ...
            'xi', 3e-3);
    elseif swap == 4 % Circle with fft and ellipse
        obstacle = 'circFFtEl';
        tm = getObst(16);
        tm = rmfield(rmfield(tm, 'xi'), 'bc');
        tm = rmfield(rmfield(tm, 'grad'), 'cur');
        tm = rmfield(rmfield(tm, 'normal'), 'derGrad');
        tm.ppw = 7;
        par = struct('obsts', [struct('par', @(t) repmat([-0.5;-2], 1, length(t)) + [0.7*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*sqrt(0.7^2+1/4)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 7, 'serpar', ...
            @(t,mo) [(-0.5*((0:mo-1) == 0) + 0.7*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ... 
            (-2*((0:mo-1) == 0) +0.5*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))]) tm], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 6e-3);
    elseif swap == 5 % Circle without fft and ellipse
        obstacle = 'circWoFFtEl';
        par = struct('obsts', [struct('par', @(t) repmat([-0.5;-2], 1, length(t)) + [0.7*cos(2*pi*t); 0.5*sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*sqrt(0.7^2+1/4)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 7, 'serpar', ...
            @(t,mo) [(-0.5*((0:mo-1) == 0) + 0.7*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ... 
            (-2*((0:mo-1) == 0) +0.5*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2))])   ...
            struct('par', @(t) [cos(2*pi*t); sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 7, 'serpar', ...
            @(t,mo) [ (2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2); ... 
            (2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)]) ], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 6e-3);
    elseif swap == 6 % Near-inclusion near convex part, near-convex and ellipse
        obstacle = 'nincNcEll';
        tmp = getObst(4, [+0.5; -1]);
        tmp = rmfield(rmfield(tmp, 'xi'), 'bc');
        tmp = rmfield(rmfield(tmp, 'grad'), 'cur');
        tmp = rmfield(rmfield(tmp, 'normal'), 'derGrad');
        tmp.ppw = 10;
        tm = getObst(3, [+0; +0.6]);
        tm = rmfield(rmfield(tm, 'xi'), 'bc');
        tm = rmfield(rmfield(tm, 'grad'), 'cur');
        tm = rmfield(rmfield(tm, 'normal'), 'derGrad');
        tm.ppw = 10;
        par = struct('obsts', [tm  tmp  struct('par', @(t) repmat([-0.4;0], 1, length(t)) + [0.15*cos(2*pi*t); 0.25*sin(2*pi*t)], 'gradnorm',...
            @(t) 2*pi*sqrt(0.15^2+0.25^2)*ones(size(t)), 'corners', [], 'dbf', 1, 'ppw', 10, 'serpar', ...
            @(t,mo) [(-0.4*((0:mo-1) == 0) + 0.15*(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)); ...
            0.25*(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2)]) ], ...
            'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
    end
    par.k = ks(ki);
    par.N = 0;
    par.r = zeros(2,length(par.obsts)); % ranges
    for obst = 1:length(par.obsts)
        par.obsts(obst).k = par.k;
        par.obsts(obst).N = par.obsts(obst).ppw*par.k;
        par.r(1,obst) = par.N+1;
        par.N = par.N + par.obsts(obst).N;
        par.r(2,obst) = par.N;
        
        par.obsts(obst).t = linspace(0,1,par.obsts(obst).N+1); % The knots of the periodic spline;
        par.obsts(obst).colltau = par.obsts(obst).t(1:par.obsts(obst).N);
        par.obsts(obst).hs = (par.obsts(obst).colltau(2) -par.obsts(obst).colltau(1) )/2;
    end
    
    if bc == 4
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
    elseif bc == 3
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
    elseif bc == 2
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
    elseif bc == 1
        par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
    end
    
    overs1 = 1;
    
    %% Computating full solution
    A1 = zeros(overs1*par.N, par.N);
    b1 = zeros(overs1*par.N,1);
    tic
    prevToc = toc;
    for i = 1:par.N
        if (toc-prevToc > printtoc)
            prevToc = toc;
            display([num2str(par.k), '=k, ' num2str( (i-1)/par.N,'%7.3f'), '=A1%, est. # sec. left for A1=' num2str(toc*(par.N-i+1)/(i-1)) ])
        end
        obsin = find((i >= par.r(1,:)) & (i <= par.r(2,:)));
        collx = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) );
        
        collxos = zeros(2,overs1);
        for os = 1:overs1
            collxos(:,os) = par.obsts(obsin).par(par.obsts(obsin).colltau(i-par.r(1,obsin)+1) + par.obsts(obsin).hs*2/overs1*(os-1));
        end
        for obst = 1:length(par.obsts)
            if obst == obsin
                for os = 1:overs1
                    A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF(i-par.r(1,obst)+1,par.obsts(obst), 1,par.obsts(obst).N, [], [], ...
                        par.obsts(obst).hs*2/overs1*(os-1) );
                    b1(overs1*i -overs1+os) = par.bc(par.k, collxos(:,os));
                end
            else
                for os = 1:overs1
                    A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                end
            end
        end
    end
    c1 = A1\b1;
    % Possibly compute and plot the correlations
    % [R, sigma,obbounds] = calcCorr(par, c1, Tcor, percDecay);
    % [R, sigma,obbounds] = calcCorr(par, c1, 0.1, 1, [10, now], A1);
    
    % figure; spectrogram(c1, round(length(c1)/16),[],[], 1/(par.obsts(1).colltau(2)-par.obsts(1).colltau(1)), 'centered'); 
    % figure; plot( [real(c1) imag(c1)]); legend('Re(c1)','Im(c1)')
    % v = validate(A1,nan*A1,par,v,idx)
    
    %% Eigenvalues big reflection matrix
    bigm = eye(par.obsts(1).N);
    for obst = 1:length(par.obsts)
        i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
        % Need brackets such that the system matrix is A(o,o) iso all previous times A(o,o) -> bad cond
        bigm = (A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\A1(...
            par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst)))*bigm;
    end
    display(['Starting EVD of matrix of size ' num2str(size(bigm))]);
    [Vb, Db, Wb] = eig(bigm);
    display(['Ended EVD of matrix of size ' num2str(size(bigm))]);
    
    % Check whether the eigenvalue is as expected; modulus also for negative angles
    [real(Db(1,1)) imag(Db(1,1)) abs(Db(1,1)) mod(angle(Db(1,1)),2*pi) mod(2*1*par.k,2*pi)]
    
    
    %% Saving 'mev' eigenvectors per obstacle
    V1 = Vb(:,1);
    mev = 2;
    allV = cell(length(par.obsts)+1,1);
    allV{1} = Vb(:,1:mev);
    for obst = 1:length(par.obsts)
        matr = zeros(par.obsts(obst).N, mev);
        for ev = 1:mev
                i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
                matr(:,ev) = A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\(A1(...
                    par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst))*allV{obst}(:,ev));
        end
        allV{obst+1} = matr;
    end
    for ev = 1:mev
        if norm(allV{1}(:,ev) - allV{length(par.obsts)+1}(:,ev)/Db(ev,ev) ) > norm(allV{1}(:,ev))*eps^(4/7) 
            error('Not an eigenvector');
        end
    end
    save(['V1k' num2str(par.k) 'obst' num2str(obstacle) '.mat'], 'V1', 'par', 'c1', 'allV');
    
    display(['ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2) )  ]);
end

