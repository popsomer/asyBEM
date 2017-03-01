% Make plots of the eigenvector, phitilde and convergence plots for two circles.

%% Loading saved eigenvector
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

load V1k9; % Loads V1l and par
obst = 1;

% Make sure tau_1 = 0 is in the middle
l = length(par.obsts(obst).colltau)/2;
signal = [V1(l+1:end); V1(1:l)];
collsignal = [(par.obsts(obst).colltau(l+1:end)-1), par.obsts(obst).colltau(1:l)];

rang = find(abs(collsignal) < 0.26);


%% Geometric attempts for the phase
zeta = transpose(sqrt( sum( (par.obsts(obst).par(collsignal) - repmat(par.obsts(3-obst).par(0), 1, ...
    length(collsignal) ) ).^2, 1) ) );

c1x = 0; c1y = 0; c2x = 0; c2y = 2; r1 = 0.5; r2 = 0.5; d = 1;
tau2 = asin( (cos(2*pi*collsignal) -4)/2./sin(2*pi*collsignal) + ...
    (-1).^(collsignal < 0)./2.*sqrt((cos(2*pi*collsignal) -4).^2./sin(2*pi*collsignal).^2 +4) )/2/pi;
xi = transpose(sqrt( sum( (par.obsts(1).par(collsignal) - par.obsts(2).par(tau2) ).^2, 1) ) );


%% Symbolic results for the phase
t1m = transpose(collsignal);
symbTay = d*ones(length(collsignal), 5);
symbTay(:,2) = d + sqrt(2)*pi^2*t1m.^2;
symbTay(:,3) = symbTay(:,2) -11/12*sqrt(2)*pi^4*t1m.^4;
symbTay(:,4) = symbTay(:,3) + 2783/2520*sqrt(2)*pi^6*t1m.^6;
symbTay(:,5) = symbTay(:,4) -358021/205632*sqrt(2)*pi^8*t1m.^8;
% Symbolic $c_i$ & $1$  & $\sqrt{2}\pi^2$ & $\frac{-11}{12}\sqrt{2}\pi^4$ & $\frac{2783\sqrt{2}\pi^6}{2520}$ & $\frac{-358021}{205632}\sqrt{2}\pi^8$ \\


%% Calculate phitilde: extract phase from angle
phitilde = zeros(size(signal));
obst = 1;

closest = find(abs(collsignal-0) == min(abs(collsignal-0)));
phitilde(closest) = 1;
for i = (closest-1):-1:1
    phitilde(i) = phitilde(i+1) +(angle(signal(i)) - angle(signal(i+1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i+1)))/2/pi);
end
for i = (closest+1):length(V1)
    phitilde(i) = phitilde(i-1) +(angle(signal(i)) - angle(signal(i-1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i-1)))/2/pi);
end

figure; plot(collsignal(rang), phitilde(rang), 'g'); hold on
plot(collsignal(rang), zeta(rang), 'b--');
plot(collsignal(rang), xi(rang), 'r:');
legend({'$\tilde{\phi}$', '$\zeta$', '$\xi$'}, 'interpreter', 'latex', 'FontSize', 23);
xlabel('\tau_1')
set(gca,'FontSize',20)


%% Plot eigenvector
figure; 

if 1
    % Plot for all tau_1
    V1plc8 = plot(collsignal, real(signal./exp(1i*par.k*symbTay(:,5))), 'g', 'LineWidth', 3);
    hold on;
    V1plc2 = plot(collsignal, real(signal./exp(1i*par.k*symbTay(:,2))), 'r--', 'LineWidth', 2);
    V1pl = plot(collsignal, real(signal), 'b:');
else
    % Plot only the region where the circles can approximately `see' each other
    V1plc8 = plot(collsignal(rang), real(signal(rang)./exp(1i*par.k*symbTay(rang,5))), 'g', 'LineWidth', 3);
    hold on;
    V1plc2 = plot(collsignal(rang), real(signal(rang)./exp(1i*par.k*symbTay(rang,2))), 'r--', 'LineWidth', 2);
    V1pl = plot(collsignal(rang), real(signal(rang)), 'b:');
end

legend([V1pl, V1plc2, V1plc8], {'Re($V_{j,1})$', 'Re$\{V_{j,1}/\exp(ik[c_0+c_2\tau^2])\}$', 'Re$\{V_{j,1}/\exp(ik\sum_{i=0}^8 c_i \tau_1^i)\}$'},...
    'interpreter', 'latex', 'FontSize', 15);

%% Plot convergence
lws = 'LineWidth';
lw = 2;
if 1
    figure; 
    loglog(abs(collsignal(rang)), (zeta(rang) -phitilde(rang))./phitilde(rang), 'b--', lws, lw);
    hold on;
    loglog(abs(collsignal(rang)), (phitilde(rang) -xi(rang))./phitilde(rang), 'r:', lws, lw);
    loglog(abs(collsignal(rang)), (phitilde(rang) -symbTay(rang,1))./phitilde(rang), 'g', lws, lw);
    loglog(abs(collsignal(rang)), abs(-symbTay(rang,2) +phitilde(rang))./phitilde(rang), 'm-.', lws, lw);
    loglog(abs(collsignal(rang)), (phitilde(rang) -symbTay(rang,3))./phitilde(rang), 'c--', lws, lw);
    loglog(abs(collsignal(rang)), abs(symbTay(rang,4) -phitilde(rang))./phitilde(rang), 'y', lws, lw);
    loglog(abs(collsignal(rang)), (phitilde(rang) -symbTay(rang,5))./phitilde(rang), 'k:', lws, lw);
else
    figure
    semilogy(collsignal, (zeta - phitilde)./(phitilde), 'b');
    hold on;
    semilogy(collsignal, abs(xi-phitilde)./phitilde, 'r');
    semilogy(collsignal, (phitilde -symbTay(:,1))./phitilde, 'g');
    semilogy(collsignal, abs(-symbTay(:,2)+phitilde)./phitilde, 'm');
    semilogy(collsignal, abs(symbTay(:,3)-phitilde)./phitilde, 'c');
    semilogy(collsignal, abs(symbTay(:,4)-phitilde)./phitilde, 'y');
    semilogy(collsignal, abs(symbTay(:,5)-phitilde)./phitilde, 'k');
end

legend({'$(\zeta -\tilde{\phi})/\tilde{\phi}$', '$(\tilde{\phi}-\xi)/\tilde{\phi}$',...
    '$(\tilde{\phi} -c_0)/\tilde{\phi}$','$|c_0+c_2\tau_1^2-\tilde{\phi}|/\tilde{\phi}$',...
    '$(\tilde{\phi} -\sum_{i=0}^4 c_i \tau_1^i)/\tilde{\phi}$','$|-\tilde{\phi}+\sum_{i=0}^6 c_i \tau_1^i|/\tilde{\phi}$',...
    '$(\tilde{\phi} -\sum_{i=0}^8 c_i \tau_1^i)/\tilde{\phi}$'}, 'interpreter', 'latex', 'FontSize', 15);
xlabel('\tau_1')
set(gca,'FontSize',13)
ylabel('Relative error')


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


%% Attempts to solve the nonlinear system of PDEs
nr = 1e3;
cl1 = linspace(0,1/2,nr);
h = cl1(2)-cl1(1);
clos = @(phi,tau2, shft) interp1(cl1-shft*h, phi, tau2);
dist = @(tau2) sqrt( (c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*tau2) ).^2 + (c1y-c2y+r2*sin(2*pi*cl1)-r2*sin(2*pi*tau2) ).^2 );

F = @(x) [ (x(1:nr)-dist(x(nr+1:end))-clos(x(1:nr), x(nr+1:end)-1/2, 0)-d ), ...
    ((clos(x(1:nr),x(nr+1:end)-1/2, 1)-clos(x(1:nr),x(nr+1:end)-1/2, 0))/h + ...
    2*pi*(r1*sin(2*pi*x(nr+1:end) ).*(c1x-c2x+r1*cos(2*pi*cl1)-r1*cos(2*pi*x(nr+1:end))) ...
    -r2*cos(2*pi*x(nr+1:end)).*(c1y-c2y+r2*sin(2*pi*cl1)+r2*sin(2*pi*x(nr+1:end))) )./dist(x(nr+1:end) ) ) ];

x0 = [sqrt( sum( (par.obsts(1).par(cl1) - repmat(par.obsts(2).par(0.75), 1, length(cl1) ) ).^2, 1) ), 0.75*ones(size(cl1))];

if 0 % Using fsolve
    opto = optimoptions('fsolve', 'Display', 'iter-detailed', 'FunctionTolerance', 1e-2, 'OptimalityTolerance', 1e-3, 'StepTolerance', 1e-12);
    [X, FVAL] = fsolve(F, x0, opto);
    
    trph = sqrt( sum( (par.obsts(obst).par(cl1) - repmat(par.obsts(2).par(0.75), 1, nr) ).^2, 1) );
    figure; plot(cl1, X(1:nr)- trph); title('diff \phi')
    figure; plot(cl1, X((nr+1):end)); title('\tau_2')
    
elseif 0 % Iterative solve
    figure; subplot(1,2,1); plot(x0(1:nr)); hold on; subplot(1,2,2); plot(x0(nr+1:end)); hold on;
    x = F(x0); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
    x = F(x); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
    x = F(x); subplot(1,2,1); plot(x(1:nr)); subplot(1,2,2); plot(x(nr+1:end));
    title('\tau_2'); legend('x0','x1', 'x2', 'x3');
    subplot(1,2,1); title('\phi(\tau_1)'); legend('x0','x1', 'x2', 'x3');
end
