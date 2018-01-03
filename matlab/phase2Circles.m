% Make plots of the eigenvector, phitilde and convergence plots for two circles.
% addpath('/dir/to/chebfun')

%% Loading saved eigenvector
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

load V1k9; % Loads V1l and par

% Make sure tau_1 = 0 is in the middle
l = length(par.obsts(1).colltau)/2;
signal = [V1(l+1:end); V1(1:l)];
collsignal = [(par.obsts(1).colltau(l+1:end)-1), par.obsts(1).colltau(1:l)];

dom = [-1/5,1/5]; % A domain smaller than the support of the mode such that the inverse of chi exists
rang = find(abs(collsignal) < 1/5);
tau = chebfun('tau',dom);
evalCh = @(f,g) f(g);


%% Geometric attempts for the phase

c1x = 0; c1y = 0; c2x = 0; c2y = 2; r1 = 0.5; r2 = 0.5; d = 1;
if 1
    zeta = transpose(sqrt( sum( (par.obsts(1).par(collsignal) - repmat(par.obsts(2).par(0), 1, ...
        length(collsignal) ) ).^2, 1) ) );
    tau2Xi = asin( (cos(2*pi*collsignal) -4)/2./sin(2*pi*collsignal) + ...
        (-1).^(collsignal < 0)./2.*sqrt((cos(2*pi*collsignal) -4).^2./sin(2*pi*collsignal).^2 +4) )/2/pi;
    xi = transpose(sqrt( sum( (par.obsts(1).par(collsignal) - par.obsts(2).par(tau2Xi) ).^2, 1) ) );
else
    % Chebfun approximation
    zeta = sqrt((evalCh(par.obsts(1).par(tau), 1) - evalCh(par.obsts(2).par(0), 1) ).^2 + ...
        (evalCh(par.obsts(1).par(tau), 2) - evalCh(par.obsts(2).par(0), 2) ).^2);
    zeta = zeta.blocks{1,1};
    
    tau2XiPts = asin( (cos(2*pi*collsignal) -4)/2./sin(2*pi*collsignal) + ...
        (-1).^(collsignal < 0)./2.*sqrt((cos(2*pi*collsignal) -4).^2./sin(2*pi*collsignal).^2 +4) )/2/pi;
    tau2XiPts(isnan(tau2XiPts)) = 0; % This should be at collsignal == 0
    tau2Xi = chebfun(transpose(tau2XiPts), 'equi');
    tau2Xi.domain = dom;
    tau2Xi.funs{1,1}.domain = tau2Xi.domain;
    xi = sqrt((evalCh(par.obsts(1).par(tau), 1) - sin(2*pi*tau2Xi)/2 ).^2 + ...
        (evalCh(par.obsts(1).par(tau), 2) - (2 -cos(2*pi*tau2Xi)/2) ).^2);
    xi = xi.blocks{1,1};
end


%% Symbolic results for the phase
ci = [d, 0, (sqrt(2)*pi^2), 0, (-11/12*sqrt(2)*pi^4), 0, (2783/2520*sqrt(2)*pi^6), 0, (-358021/205632*sqrt(2)*pi^8)];
% Symbolic $c_i$ & $1$  & $\sqrt{2}\pi^2$ & $\frac{-11}{12}\sqrt{2}\pi^4$ & $\frac{2783\sqrt{2}\pi^6}{2520}$ & $\frac{-358021}{205632}\sqrt{2}\pi^8$ \\
symbTay0 = ci(1) +0*tau;
symbTay2 = ci(1) + ci(3)*tau.^2;
symbTay4 = symbTay2 +ci(5)*tau.^4;
symbTay6 = symbTay4 + ci(7)*tau.^6;
symbTay8 = symbTay6 + ci(9)*tau.^8;

% tau2 for ray that ends up at tau1=a1^2*tau1 + ... after the next reflection so that it can reach tau1 =0 =tau2 after 
% an infinite number of reflections for |a1|<1, infTay is the series of the distance to the next SP
infTay0 = d +0*tau;
infTay2 = d + 1/2*pi^2*(3*(2*sqrt(2) - 3)^2 + 4*sqrt(2) - 3)*tau.^2;
infTay4 = infTay2 + ( -1/24*pi^4*(39*(2*sqrt(2) - 3)^4 + 52*(2*sqrt(2) - 3)^3 + 42*(2*sqrt(2) - 3)^2 + 104*sqrt(2) - 117)...
    + 7*pi^2*(3*pi^2*(17*sqrt(2)- 24)*(2*sqrt(2) - 3) + pi^2*(17*sqrt(2) - 24)) )*tau.^4;

% Test whether the symbolic results are the ones computed by the general algorithm
[taus, c, a, ft] = seriesPhasePerOrbit(getObst(5), 8);
[norm(c(:,1)), norm(c(:,2)-sqrt(2)*pi^2), norm(c(:,3)), norm(c(:,4) +11/12*sqrt(2)*pi^4), norm(c(:,5)), ...
    norm(c(:,6) - 2783/2520*sqrt(2)*pi^6), norm(c(:,7)), norm(c(:,8) + 358021/205632*sqrt(2)*pi^8)]
'= abs error on symbolic c, abs error on symbolic a ='
[norm(a(:,1,1) - (3-2*sqrt(2))), norm(a(:,1,2)), norm(a(:,1,3) +7*(17*sqrt(2)-24)*pi^2), norm(a(:,1,4)), ...
    norm(a(:,1,5) +(1205811*sqrt(2) -1705312)*pi^4/84), norm(a(:,1,6)), norm(a(:,1,7,:)  + pi^6/128520*(289615597399*sqrt(2) -409578202752))]


%% Calculate phitilde: extract phase from angle
phitilde = zeros(size(signal));

closest = find(abs(collsignal-0) == min(abs(collsignal-0)));
phitilde(closest) = 1;
for i = (closest-1):-1:1
    phitilde(i) = phitilde(i+1) +(angle(signal(i)) - angle(signal(i+1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i+1)))/2/pi);
end
for i = (closest+1):length(V1)
    phitilde(i) = phitilde(i-1) +(angle(signal(i)) - angle(signal(i-1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i-1)))/2/pi);
end


%% Plot eigenvector
figure;
V1plc8 = plot(collsignal, real(transpose(signal)./exp(1i*par.k*symbTay8(collsignal) )), 'g', 'LineWidth', 3);
hold on;
V1plc2 = plot(collsignal, real(transpose(signal)./exp(1i*par.k*symbTay2(collsignal))), 'r--', 'LineWidth', 2);
V1pl = plot(collsignal, real(transpose(signal)), 'b:');
legend([V1pl, V1plc2, V1plc8], {'Re$\{\tilde{V}(\tau)\}$', 'Re$\{\tilde{V}(\tau)/\exp(ik[c_0+c_2\tau^2])\}$', ...
    'Re$\{\tilde{V}(\tau)/\exp(ik\sum_{i=0}^8 c_i \tau^i)\}$'},...
    'interpreter', 'latex', 'FontSize', 20); %15);
xlabel('\tau');
ylabel('Mode');
xticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]);
set(gca, 'FontSize', 20);

%% Compute the extrapolated reference result
if 0
    rg = round(length(collsignal)/4);
    figure; plot(collsignal(rg:end-rg), phitilde(rg:end-rg)); title('phitilde'); % -> D268 Basic fitting to get data below
    data = [9.725, 12.85, 13.72, 13.911, 13.951]; figure; plot(data); hold on; % For c2
    % data = [-eps, -58.312, -99.94, -117.2, -122.69]; figure; plot(data); hold on; % For c4
    % data = [488.41, 966.1, 1231.7]; figure; plot(data); hold on; % For c6
    
    sgn = sign(data(1));
    if norm(sgn-sign(data)) ~= 0, error('Signs not equal'); end
    
    abc = lsqnonlin(@(abc) abc(3) -sgn*exp(abc(1)-(1:length(data))*abc(2)) -data, [1,1,14]) % Last entry is mu
    
    plot(abc(3) -sgn* exp(abc(1)-(1:length(data))*abc(2)));
end


%% Solve the nonlinear equation using chebfun
chi = (3-2*sqrt(2))*tau; % Initial guess
F = @(g, tau1) 4*pi*g(tau1) + atan( (0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*tau1)   )./(1+2*0.5 -0.5*cos(2*pi*g(tau1)) -0.5*cos(2*pi*tau1) ) ) + ...
    atan((0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*g(g(tau1))) )./(1+2*0.5 -0.5*cos(2*pi*g(tau1)) -0.5*cos(2*pi*g(g(tau1)) )) );

dF = @(g, tau1, dg) 4*pi + 1./(1+ ((0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*tau1) )./( 1+2*0.5 -0.5*cos(2*pi*g(tau1)) -0.5*cos(2*pi*tau1) ) ).^2).*...
    (0.5*2*pi.*cos(2*pi*g(tau1)).*( 1+2*0.5 -0.5*cos(2*pi*g(tau1)) -0.5*cos(2*pi*tau1) )  ...
    - (0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*tau1) ).*(-0.5)*(-2*pi).*sin(2*pi*g(tau1))  )./( 1+2*0.5 -0.5*cos(2*pi*g(tau1)) -0.5*cos(2*pi*tau1) ).^2 + ...
    1./(1 +((0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*g(g(tau1))) )./(1+2*0.5 -0.5*cos(2*pi*g(tau1))  -0.5*cos(2*pi*g(g(tau1)) )) ).^2).*...
    ( (0.5*2*pi*cos(2*pi*g(tau1)) -0.5*2*pi*cos(2*pi*g(g(tau1))).*dg(g) ).*(1+2*0.5 -0.5*cos(2*pi*g(tau1))  -0.5*cos(2*pi*g(g(tau1)) )) -...
    (0.5*sin(2*pi*g(tau1)) -0.5*sin(2*pi*g(g(tau1))) ).*( 0.5*2*pi*sin(2*pi*g(tau1)) +0.5*2*pi*sin(2*pi*g(g(tau1)) ).*dg(g)  )...
    )./(1+2*0.5 -0.5*cos(2*pi*g(tau1))  -0.5*cos(2*pi*g(g(tau1)) )).^2;

nbIt = 6;
normFchi = zeros(nbIt+1,1);
normFchi(1) = norm(F(chi,tau));
for r = 1:nbIt
    chi = chi -F(chi,tau)./dF(chi,tau,diff(chi));
    normFchi(r+1) = norm(F(chi,tau));
end
normFchi % Machine precision after 6 iterations

endI = 9;
as =  [0, (3-2*sqrt(2)), 0, 7*pi^2*(24-17*sqrt(2)), 0, -1/84*pi^4*(1205811*sqrt(2) - 1705312), 0, ...
    -1/128520*pi^6*(289615597399*sqrt(2) - 409578202752), zeros(1, endI-7)];
converAi = zeros(endI+1,4);
for i = 0:endI
    converAi(i+1,:) = [i, abs((evalCh(diff(chi,i), 0)/gamma(i+1)- as(i+1))/as(i+1)), evalCh(diff(chi,i), 0), as(i+1)];
end
converAi % Compare chi to its series expansion


%% Compute the phase
tau2 = chebfun('tau2',dom);
dist = chebfun2( @(tau,tau2) sqrt( (0.5*sin(2*pi*tau2) -0.5*sin(2*pi*tau) ).^2 + (1+2*0.5 -0.5*cos(2*pi*tau2) -0.5*cos(2*pi*tau) ).^2 ), ...
    [-1, 1, -1, 1]/5);
ddist = diffy(dist);

phi = dist(tau,chi) - 1;
prev = tau;
chiRefl = chi;

nbRefl = 10;
cfs = zeros(length(ci),nbRefl+1);
for i = 0:length(ci)-1
    cfs(i+1,1) = evalCh(diff(phi,i),0)/gamma(i+1);
end
conver = zeros(1,nbRefl+1);
conver(1) = norm(phi);
for r = 1:nbRefl
    prev = chiRefl;
    chiRefl = chi(chiRefl);
    phi = phi + dist(prev,chiRefl) -1;
    for i = 0:length(ci)-1
        cfs(i+1,r+1) = evalCh(diff(phi,i),0)/gamma(i+1);
    end
    conver(r+1) = norm( dist(prev,chiRefl) -1);
end
phi = phi + ci(1);
conver % Show the difference in successive approximations
converCi = (cfs - repmat(transpose(ci), 1, nbRefl+1))./repmat(transpose(ci), 1, nbRefl+1) % Show the convergence to the series expansion coefficients

phiInt = evalCh(cumsum(-ddist(tau,chi).*diff(chi)), inv(chi));
phiInt = phiInt - phiInt(0) + ci(1);
phiRestr = restrict(phi, phiInt.domain);


%% Plot phi
figure; plot(phi, 'g'); hold on
plot(collsignal(rang), zeta(rang), 'b--');
plot(collsignal(rang), xi(rang), 'r:');

% Issue 2061 on github/chebfun
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'g');
h(2) = plot(NaN,NaN,'b--');
h(3) = plot(NaN,NaN,'r:');
legend(h,{'$\tilde{\phi}$', '$\zeta$', '$\xi$'}, 'interpreter', 'latex', 'FontSize', 23);
xlabel('\tau_1')
set(gca,'FontSize',20)


%% Plot the convergence
lws = 'LineWidth';
lw = 2;
chr = [1e-4, max(dom)];

figure;
loglog(abs(collsignal(rang)), (zeta(rang) -phitilde(rang))./phitilde(rang), 'b', lws, lw);
hold on;
loglog(abs(collsignal(rang)), (phitilde(rang) -xi(rang))./phitilde(rang), 'r:', lws, lw);
loglog((phi -symbTay0)./phi, chr, 'g', lws, lw);
loglog((symbTay2 -phi)./phi, chr, 'm-.', lws, lw);
loglog((phi -symbTay4)./phi, chr, 'c--', lws, lw);
loglog(abs(symbTay6 -phi)./phi, chr, 'y', lws, lw);
loglog(abs(phi -symbTay8)./phi, chr, 'k:', lws, lw);

loglog(abs(collsignal(rang)), abs(transpose(phitilde(rang)) -phi(collsignal(rang)))./phi(collsignal(rang)), 'k', lws, lw);
loglog(abs(phiRestr - phiInt)./phiRestr, 'g--', lws, lw);

h = zeros(9, 1);
h(1) = plot(NaN,NaN,'g');
h(2) = plot(NaN,NaN,'b');
h(3) = plot(NaN,NaN,'r:');
h(4) = plot(NaN,NaN,'k');
h(5) = plot(NaN,NaN,'m-.');
h(6) = plot(NaN,NaN,'c--');
h(7) = plot(NaN,NaN,'y');
h(8) = plot(NaN,NaN,'k:');
h(9) = plot(NaN,NaN,'g--');
legend(h, {'$(\phi -c_0)/\phi$', '$(\zeta -\phi)/\phi$', '$(\phi-\xi)/\phi$',...
    '$|\phi-\tilde{\phi}|/\phi$', '$(c_0+c_2\tau_1^2-\phi)/\phi$',...
    '$(\phi -\sum_{i=0}^4 c_i \tau_1^i)/\phi$','$|-\phi+\sum_{i=0}^6 c_i \tau_1^i|/\phi$',...
    '$|\phi -\sum_{i=0}^8 c_i \tau_1^i|/\phi$',...
    '$|\phi + \int \delta''(\tau_1,\chi) \chi'' d\tau_1|/\phi$'}, 'interpreter', 'latex', ...
    'FontSize', 20, 'Location', 'best'); %'southeastoutside');
%     '$|\phi + \int \delta''(\tau_1,\chi) \chi'' d\tau_1|/\tilde{\phi}$'}, 'interpreter', 'latex', 'FontSize', 15);
ylim([eps,1]);
xlim(chr);
xlabel('\tau');
set(gca,'FontSize', 20); %13);
ylabel('Relative error');
