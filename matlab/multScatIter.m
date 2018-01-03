% Compute the iterations of the reflections in a multiple scattering configuration.

%% Initialising
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

ks = 2^7;
obstacle = 5; % Two circles

bc = 5;

printtoc = 10;
kl = length(ks);
nbOb = 1;

avm = 100; % Number of random taus to average BC over
v = struct('avm', avm, 'taus', rand(avm,1), 'errBCavm', zeros(nbOb*kl,2), 'errSol', zeros(nbOb*kl,2),  ...
    'errInt', zeros(nbOb*kl,2), 'timeA', zeros(nbOb*kl,4), 'ks', ks); %, 'field', zeros(70) );

%% Computations
start = now;
ki = 1;
par = getObst(obstacle); % Reset par

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

if bc == 5
    par.bc = @(k,x) -1*besselh(0,1, k*sqrt( (x(1,:)' - 0.5).^2 + (x(2,:)' -0.5).^2) ); % Pt source
elseif bc == 4
    par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2-bcsh);sin(-pi/2-bcsh)]);
elseif bc == 3
    par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/2);sin(-pi/2)]);
elseif bc == 2
    par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(-pi/4);sin(-pi/4)]);
elseif bc == 1
    par.bc = @(k,x)-1*exp(1i*k*(x')*[cos(0);sin(0)]);
end


%% Computating full solution
overs1 = 1; % if > 1 then oversampling
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

%% Multiple scattering iterations by extraction of submatrices
nbRefl = 3;
cs = cell(length(par.obsts), nbRefl+1);

legs = cell(nbRefl+1,1);
legs{1} = 'From incident wave';
marks = {'b:', 'r-.', 'k', 'm'};
endo = length(par.obsts);

% Compute the density associated to the incident wave
for obst = 1:endo
    i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
    cs{obst,1} = A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\b1(par.r(1,i):par.r(2,i));
    figure(i); hold on;
    toShft = par.obsts(i).colltau > 0.5;
    plot([(par.obsts(i).colltau(toShft)-1), par.obsts(i).colltau(~toShft)], ...
        [real(cs{obst,1}(toShft)); real(cs{obst,1}(~toShft))], marks{1}, 'LineWidth', 2);
    xlabel(['\tau_' num2str(i)]); ylabel(['Re [v_' num2str(i) '(\tau_' num2str(i) ')]']);
end

% Compute consecutive reflections
for refl = 1:nbRefl
    for obst=1:endo
        i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
        legs{refl+1} = ['Reflection ' num2str(refl)];
        cs{obst,refl+1} = (A1(par.r(1,i):par.r(2,i), par.r(1,i):par.r(2,i))\A1(...
            par.r(1,i):par.r(2,i),par.r(1,obst):par.r(2,obst)))*cs{obst,refl}; %cs{i,refl}
        figure(i);
        toShft = par.obsts(i).colltau > 0.5;
        plot([(par.obsts(i).colltau(toShft)-1), par.obsts(i).colltau(~toShft)], [real(cs{obst, refl+1}(toShft)); ...
            real(cs{obst, refl+1}(~toShft))], marks{refl+1}, 'LineWidth', 2+refl/nbRefl);
    end
end
for obst = 1:endo
    i = obst +1 -(obst == length(par.obsts))*length(par.obsts);
    figure(i);
    legend(legs, 'location', 'best');
    set(gca, 'FontSize', 20);
    xticks([-0.5, -0.4, -0.25, -0.1, 0, 0.1, 0.25, 0.4, 0.5]);
    if bc == 5
        axis([-1/4, 1/4, -30, 30])
        xticks([-0.25, -0.2, -0.125, -0.05, 0, 0.05, 0.125, 0.2, 0.25]);
    end
end
hold off
