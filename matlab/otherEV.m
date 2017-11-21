% Compute the influence of the number of points per wavelength, the wavenumber and the geometry on all (significant) 
% EigenValues and V_{j,i}, the i-th eigenvector of M which represents a full cycle of reflections between the two circles.

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^7; %2^8;
d = 1; %2; %4; %Dist between circles
% Could also adjust radius r but scales with d and k
ppw = 5; %10; %20;

obstacle = 5; % Two circles
printtoc = 10;
kl = length(ks);
nEV = 30; %6;


%% Computations
start = now;

for ki = 1:kl % Or vary d or ppw in this loop
    par = getObst(obstacle);
    
    par.obsts(2).par = @(t) repmat([0; 1+d], 1, length(t)) + [sin(2*pi*t); -cos(2*pi*t)]/2;
    par.obsts(2).serpar = @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ...
            ((2+2*d)*((0:mo-1) == 0) -(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2;
%     par = struct('obsts', [struct('par', @(t) [sin(2*pi*t); cos(2*pi*t)]/2, 'gradnorm',... 
% 			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
%             'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ... 
%             (2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2)]/2) ...
% 			struct('par', @(t) repmat([0; 2], 1, length(t)) + [sin(2*pi*t); -cos(2*pi*t)]/2, 'gradnorm',...
% 			@(t) pi*ones(size(t)), 'normal', @(t) [sin(2*pi*t); -cos(2*pi*t)], 'corners', [], 'dbf', 1, 'ppw', 10, ...
%             'serpar', @(t,mo) [(2*pi).^(0:mo-1)./factorial(0:mo-1).*sin(2*pi*t + (0:mo-1)*pi/2); ...
%             (4*((0:mo-1) == 0) -(2*pi).^(0:mo-1)./factorial(0:mo-1).*cos(2*pi*t + (0:mo-1)*pi/2))]/2)], ...
%             'bc', @(k,x) -1*exp(1i*k*(x')*[cos(0); sin(0)]), 'xi', 3e-3);
        
    par.k = ks(ki);
    par.N = 0;
    par.r = zeros(2,length(par.obsts)); % ranges
    for obst = 1:length(par.obsts)
        par.obsts(obst).ppw = ppw;
        
        par.obsts(obst).k = par.k;
        par.obsts(obst).N = par.obsts(obst).ppw*par.k;
        par.r(1,obst) = par.N+1;
        par.N = par.N + par.obsts(obst).N;
        par.r(2,obst) = par.N;
        
        par.obsts(obst).t = linspace(0,1,par.obsts(obst).N+1); % The knots of the periodic spline;
        par.obsts(obst).colltau = par.obsts(obst).t(1:par.obsts(obst).N);
        par.obsts(obst).hs = (par.obsts(obst).colltau(2) -par.obsts(obst).colltau(1) )/2;
    end
    
    overs1 = 1; % possible oversampling?
    
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
%                     b1(overs1*i -overs1+os) = par.bc(par.k, collxos(:,os));
                end
            else
                for os = 1:overs1
                    A1(overs1*i-overs1+os,par.r(1,obst):par.r(2,obst)) = collRowQBF('error',par.obsts(obst),1,par.obsts(obst).N, [], collxos(:,os) );
                end
            end
        end
    end
    
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
    
    
    display(['ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2) )  ]);
    
    figure; semilogy(abs(diag(Db))); ylabel('|D_{i,i}|'); xlabel('i');
    phasEV = [mod(2*d*par.k,2*pi); angle(diag(Db(1:nEV,1:nEV)))]
    
    params = [ks, d, ppw, nEV]
    for ev = 1:nEV
%         signal = Vb(:,ev);
%         collsignal = par.obsts(ev).colltau;
        signal = [Vb((par.N/4+1):par.N/2,ev); Vb(1:par.N/4,ev)];
        collsignal = [par.obsts(1).colltau((par.N/4+1):par.N/2)-1, par.obsts(1).colltau(1:par.N/4)];
        phitilde = zeros(size(signal));
        closest = find(abs(collsignal -0) == min(abs(collsignal-0))); % 0 is tau_1^*
        phitilde(closest) = angle(signal(closest))./par.k;
        
        multipl = 0;
        for i = (closest-1):-1:1
            phitilde(i) = angle(signal(i))./par.k + multipl*2*pi/par.k;
            if angle(signal(i)) > (angle(signal(i+1)) +pi)
                multipl = multipl -1;
                phitilde(i) = phitilde(i) - 2*pi/par.k;
            elseif angle(signal(i)) < (angle(signal(i+1)) -pi)
                multipl = multipl +1;
                phitilde(i) = phitilde(i) + 2*pi/par.k;
            end
        end
        multipl = 0;
        for i = (closest+1):length(signal)
            phitilde(i) = angle(signal(i))./par.k + multipl*2*pi/par.k;
            if angle(signal(i)) > (angle(signal(i-1)) +pi)
                multipl = multipl -1;
                phitilde(i) = phitilde(i) - 2*pi/par.k;
            elseif angle(signal(i)) < (angle(signal(i-1)) -pi)
                multipl = multipl +1;
                phitilde(i) = phitilde(i) + 2*pi/par.k;
            end
        end
        phitilde = sum(d) +(phitilde - phitilde(closest) );
        
        factor = max(abs(signal))/max(abs(phitilde));
        figure; plot(collsignal,[real(signal), phitilde*factor, abs(signal)]);
        xlabel(['\tau_' num2str(ev)]);
        legend(['Re EV ' num2str(ev)], ['phase * ' num2str(factor)], 'abs');
    end

%     save(['V1k' num2str(par.k) 'obst' num2str(obstacle) '.mat'], 'V1', 'par', 'c1', 'allV');
    
end

