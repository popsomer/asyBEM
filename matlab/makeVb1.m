% Compute the V_{j,1}, the first eigenvector of M which represents a full cycle of reflections between the two circles.

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% ks = 2^9; 
ks = 2^7;
% obsts = 5;
% obstacle = 5; % Two circles
% obstacle = 12; % Three circles
obstacle = 13; % Ellipse and near-convex
if 0 % Three circles switched
    obstacle = 12; par = getObst(obstacle);
    tmppar = par.obsts(2);
    par.obsts(2) = par.obsts(3);
    par.obsts(3) = tmppar;
end

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
    par = getObst(obstacle); % Reset par
%     if obstacle == 11
%         par = getObst(obstacle, rx); % Reset par
%     else
%         par = getObst(obstacle); % Reset par
%         rx = 0.5;
%     end
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
    % title(['Spectrogram ' ])%num2str(bcsh)])
    % figure; plot( [real(c1) imag(c1)]); legend('Re(c1)','Im(c1)')
    % v = validate(A1,nan*A1,par,v,idx)
    
    display('aspodijf');
    %% Eigenvalues big reflection matrix
    bigm = eye(par.obsts(1).N);
    if (length(par.obsts) == 2) && 0
        % Wrong formula which is correct due to symmetry:
        A11 = A1(1:size(A1,1)/2, 1:size(A1,2)/2);
        A21 = A1(1:size(A1,1)/2, (1+size(A1,2)/2):end); % Action of density of obst 2 on location of obst 1
        A12 = A1((size(A1,1)/2+1):end, 1:size(A1,2)/2);
        A22 = A1((size(A1,1)/2+1):end, (1+size(A1,2)/2):end);
        bigm = A11\A21*(A22\A12);
    elseif 0
        for obst = 1:length(par.obsts)
           other = obst-1;
           if other == 0
               other = length(par.obsts);
           end
%            bigm = bigm*A1(par.r(1,obst):par.r(2,obst),par.r(1,obst):par.r(2,obst))\A1(...
%                par.r(1,obst):par.r(2,obst),par.r(1,other):par.r(2,other));
           % Need brackets such that the system matrix is A(o,o) iso all previous times A(o,o)
           bigm = bigm*(A1(par.r(1,obst):par.r(2,obst),par.r(1,obst):par.r(2,obst))\A1(...
               par.r(1,obst):par.r(2,obst),par.r(1,other):par.r(2,other)));
%            bigm = bigm*A1(par.r(1,obst):par.r(2,obst), par.r(1,obst):par.r(2,obst))\A1(...
%                par.r(1,other):par.r(2,other), par.r(1,obst):par.r(2,obst));
        end
    else
        for obst = 1:length(par.obsts)
           i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
           % Need brackets such that the system matrix is A(o,o) iso all previous times A(o,o) -> bad cond
%            bigm = bigm*(A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\A1(...
           bigm = (A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\A1(...
               par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst)))*bigm;
        end
    end
    display(['Starting EVD of matrix of size ' num2str(size(bigm))]);
    [Vb, Db, Wb] = eig(bigm);
    display(['Ended EVD of matrix of size ' num2str(size(bigm))]);
    
    % Check whether the eigenvalue is as expected; modulus also for negative angles
    [real(Db(1,1)) imag(Db(1,1)) abs(Db(1,1)) mod(angle(Db(1,1)),2*pi) mod(2*1*par.k,2*pi)]
    
    
    %% Ray tracing by partitioned block
    if 0
        c1true1 = c1(1:length(c1)/2);
        c1true2 = c1(length(c1)/2+1:end);
        
        nbrefl = 1; %8;
        bsl1 = zeros(size(A1,1)/2, nbrefl);
        bsl2 = zeros(size(A1,1)/2, nbrefl);
        csl1 = zeros(size(A1,2)/2, nbrefl);
        csl2 = zeros(size(A1,2)/2, nbrefl);
        
        bsl1(:,1) = b1(1:size(A1,1)/2);
        bsl2(:,1) = b1((size(A1,1)/2+1):end);
        csl1(:,1) = A11\bsl1(:,1);
        csl2(:,1) = A22\bsl2(:,1);
        
        resid = cell(2,nbrefl);
        resid{1,1} = c1true1 - sum(csl1,2);
        figure; plot(par.obsts(1).colltau, real(csl1(:,1)) ); 
        title('csl1 refl'); hold on;
        
        for refl = 2:nbrefl
            bsl1(:,refl) = -A21*csl2(:,refl-1);
            bsl2(:,refl) = -A12*csl1(:,refl-1);
            csl1(:,refl) = A11\bsl1(:,refl);
            csl2(:,refl) = A22\bsl2(:,refl);
            
            resid{1,refl} = c1true1 - sum(csl1,2);
            plot(par.obsts(1).colltau, real(csl1(:,refl)));
        end
        legend(num2str((1:nbrefl)'));
    end
%     v = validate(A1,A1, par, struct('field', zeros(40)));
%     v = validate(A1,A1, par, v,1);


    V1 = Vb(:,1);
    mev = 2;
    allV = cell(length(par.obsts)+1,1);
    allV{1} = Vb(:,1:mev);
    for obst = 1:length(par.obsts)
        matr = zeros(par.obsts(obst).N, mev);
        for ev = 1:mev
%             if obst == 1
%                 matr(:,ev) = Vb(:,ev);
%             else
                i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
%                 matr(:,ev) = (A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\A1(...
%                     par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst)))*allV{obst}(:,ev);
                matr(:,ev) = A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\(A1(...
                    par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst))*allV{obst}(:,ev));
%                     par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst)))*allV{obst-1}(:,ev);
%             end
        end
        allV{obst+1} = matr;
%         allV{obst} = matr;
    end
    for ev = 1:mev
%         if norm(allV{1}(:,ev) - allV{length(par.obsts)+1}(:,ev)*norm(allV{1}(:,ev))/norm(allV{length(par.obsts)+1}(:,ev))) ...
        if norm(allV{1}(:,ev) - allV{length(par.obsts)+1}(:,ev)/Db(ev,ev) ) ...
                > norm(allV{1}(:,ev))*eps^(4/7) %*eps^(3/4)
            error('Not an eigenvector');
        end
    end
%     save('V1k9.mat', 'V1', 'par', 'c1');
%     save(['V1k' num2str(par.k) 'obst' num2str(obstacle) '.mat'], 'V1', 'par', 'c1');
%     save(['V1k' num2str(par.k) 'obst' num2str(obstacle) '.mat'], 'V1', 'par', 'c1', 'allV');
    save(['V1k' num2str(par.k) 'obst' num2str(obstacle) 'swap.mat'], 'V1', 'par', 'c1', 'allV');
    display(['ki = ' num2str(ki) ', now is ' datestr(now) ', expected end ' datestr(start + ...
        (now-start)*sum(ks.^2)/sum(ks(1:ki).^2) )  ]);
end

