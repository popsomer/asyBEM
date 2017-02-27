% Make plots of the eigenvector, phitilde and convergence plots for two circles.

%% Loading saved eigenvector
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

load V1k9; % Loads V1l and par
% load V1k7; % Loads V1l and par
obst = 1;

l = length(par.obsts(obst).colltau)/2;
signal = [V1(l+1:end); V1(1:l)];
collsignal = [(par.obsts(obst).colltau(l+1:end)-1), par.obsts(obst).colltau(1:l)];


%% Geometric attempts for the phase
zeta = transpose(sqrt( sum( (par.obsts(obst).par(collsignal) - repmat(par.obsts(3-obst).par(0), 1, ...
    length(collsignal) ) ).^2, 1) ) );

% c1x = -0.5; c1y = -0.5; c2x = -0.5; c2y = 1.5; r1 = rx; r2 = 0.5; d = 1;
c1x = 0; c1y = 0; c2x = 0; c2y = 2; r1 = 0.5; r2 = 0.5; d = 1;
% tau2 is wrong because for old parametrization
% tau2 = (c1y-c2y+r2*sin(2*pi*par.obsts(1).colltau) )./sqrt((c1y-c2y)^2+2*r1*sin(2*pi*par.obsts(1).colltau)*(c1y-c2y) + r1^2);
% tau2 = (-1).^(abs(par.obsts(1).colltau-0.5) < 0.25).*asin(tau2)/2/pi+1 -0.5*(abs(par.obsts(1).colltau-0.5) < 0.25);
% xi = transpose(sqrt( (c1x-c2x+r1*cos(2*pi*par.obsts(1).colltau)-r1*cos(2*pi*tau2) ).^2 + ...
%     (c1y-c2y+r2*sin(2*pi*par.obsts(1).colltau)-r2*sin(2*pi*tau2) ).^2) );
tau2 = (c1y-c2y+r2*sin(2*pi*collsignal) )./sqrt((c1y-c2y)^2+2*r1*sin(2*pi*collsignal)*(c1y-c2y) + r1^2);
tau2 = (-1).^(abs(collsignal-0.5) < 0.25).*asin(tau2)/2/pi+1 -0.5*(abs(collsignal -0.5) < 0.25);
xi = transpose(sqrt( (c1x-c2x+r1*cos(2*pi*collsignal)-r1*cos(2*pi*tau2) ).^2 + ...
    (c1y-c2y+r2*sin(2*pi*collsignal)-r2*sin(2*pi*tau2) ).^2) );


%% Symbolic results for the phase
% t1m = transpose(collsignal(1:l)-1/4);
% symbTay = d*ones(l, 5);
t1m = transpose(collsignal);
symbTay = d*ones(length(collsignal), 5);
symbTay(:,2) = d + sqrt(2)*pi^2*t1m.^2;
symbTay(:,3) = symbTay(:,2) -11/12*sqrt(2)*pi^4*t1m.^4;
symbTay(:,4) = symbTay(:,3) + 2783/2520/sqrt(2)*pi^6*t1m.^6;
symbTay(:,5) = symbTay(:,4) -358021/205632*sqrt(2)*pi^8*t1m.^8;
% Symbolic $c_i$ & $1$  & $\sqrt{2}\pi^2$ & $\frac{-11}{12}\sqrt{2}\pi^4$ & $\frac{2783\pi^6}{2520\sqrt{2}}$ & $\frac{-358021}{205632}\sqrt{2}\pi^8$ \\


%% Calculate phitilde: extract phase from angle
% phitilde = zeros(size(V1));
phitilde = zeros(size(signal));
obst = 1;

% closest = find(abs(collsignal-0.25) == min(abs(collsignal-0.25)));
closest = find(abs(collsignal-0) == min(abs(collsignal-0)));
phitilde(closest) = 1;
for i = (closest-1):-1:1
%     phitilde(i) = phitilde(i+1) +(angle(V1(i)) - angle(V1(i+1)))/par.k -2*pi/par.k*floor( angle(V1(i)) -angle(V1(i+1))/2/pi);
    phitilde(i) = phitilde(i+1) +(angle(signal(i)) - angle(signal(i+1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i+1)))/2/pi);
end
for i = (closest+1):length(V1)
%     phitilde(i) = phitilde(i-1) +(angle(V1(i)) - angle(V1(i-1)))/par.k -2*pi/par.k*floor( angle(V1(i)) -angle(V1(i-1))/2/pi);
    phitilde(i) = phitilde(i-1) +(angle(signal(i)) - angle(signal(i-1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i-1)))/2/pi);
end

% figure; plot(collsignal(1:l), phitilde(1:l), 'g'); hold on
% plot(collsignal(1:l), zeta(1:l), 'b');
% plot(collsignal(1:l), xi(1:l), 'r');
figure; plot(collsignal, phitilde, 'g'); hold on
plot(collsignal, zeta, 'b');
plot(collsignal, xi, 'r');
legend({'$\tilde{\phi}$', '$\zeta$', '$\xi$'}, 'interpreter', 'latex', 'FontSize', 15);
xlabel('\tau_1')


%% Plot eigenvector
figure; 
% plot(collsignal, [real(V1), imag(V1), real(V1./exp(1i*par.k*symbTay(:,2))), imag(V1./exp(1i*par.k*symbTay(:,5)))]); 
plot(collsignal, [imag(signal), real(signal./exp(1i*par.k*symbTay(:,2))), imag(signal./exp(1i*par.k*symbTay(:,5)))]); 
legend({'Im $V_{j,1}$', 'Re $V_{j,1}/\exp(ik[c_0+c_2\tau^2])$', 'Im $V_{j,1}/\exp(ik\sum_{i=0}^8 c_i \tau_1^i)$'},...
    'interpreter', 'latex', 'FontSize', 15);


%% Plot convergence
% figure; plot(collsignal(1:l)-1/4, signal(1:l), 'g'); hold on
% plot(collsignal(1:l)-1/4, (tryph(1:l)-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'b');
% plot(collsignal(1:l)-1/4, (transpose(phi(1:l))-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'r');
% plot(collsignal(1:l)-1/4, (symbTay-ph(1:l)/ks(ki))./(ph(1:l)/ks(ki)), 'k');
figure; 
loglog(abs(collsignal), (phitilde-zeta)./(phitilde), 'b');
hold on;
loglog(abs(collsignal), abs(xi-phitilde)./phitilde, 'r');
loglog(abs(collsignal), (phitilde -symbTay(:,1))./phitilde, 'g');
loglog(abs(collsignal), abs(-symbTay(:,2)+phitilde)./phitilde, 'm');
loglog(abs(collsignal), abs(symbTay(:,3)-phitilde)./phitilde, 'c');
loglog(abs(collsignal), abs(symbTay(:,4)-phitilde)./phitilde, 'y');
loglog(abs(collsignal), abs(symbTay(:,5)-phitilde)./phitilde, 'k');

legend({'$(\tilde{\phi}-zeta)/\tilde{\phi}$', '$| \xi-\tilde{\phi}|/\tilde{\phi}$',...
    '$(\tilde{\phi}-c_0)/\tilde{\phi}$','$|c_0+c_2\tau_1^2-\tilde{\phi}|/\tilde{\phi}$',...
    '$|-\tilde{\phi}+\sum_{i=0}^4 c_i \tau_1^i|/\tilde{\phi}$','$|-\tilde{\phi}+\sum_{i=0}^6 c_i \tau_1^i|/\tilde{\phi}$',...
    '$|-\tilde{\phi}+\sum_{i=0}^8 c_i \tau_1^i|/\tilde{\phi}$'}, 'interpreter', 'latex', 'FontSize', 15);
xlabel('\tau_1')


%% Compute the extrapolated reference result
rg = round(length(collsignal)/4);
figure; plot(collsignal(rg:end-rg), phitilde(rg:end-rg)); title('phitilde'); % -> D268 Basic fitting to get data below
data = [9.725, 12.85, 13.72, 13.911, 13.951]; %figure; plot(data); hold on; % For c2
% data = [-eps, -58.312, -99.94, -117.2, -122.69]; %figure; plot(data); hold on; % For c4
% data = [488.41, 966.1, 1231.7]; % figure; plot(data); hold on; % For c6

sgn = sign(data(1));
if norm(sgn-sign(data)) ~= 0, error('Signs not equal'); end

abc = lsqnonlin(@(abc) abc(3) -sgn*exp(abc(1)-(1:length(data))*abc(2)) -data, [1,1,14]) % Last entry is mu

%plot(abc(3) -sgn* exp(abc(1)-(1:length(data))*abc(2)));



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
