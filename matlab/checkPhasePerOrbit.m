% Make convergence plots for the phase of the periodic orbit of general obstacles.
% addpath('/dir/to/chebfun')

%% Loading saved eigenvector
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');
 
% load V1k128obst12; % Loads V1l and par
load V1k128obst13; % Ellipse and nearconvex
if 0
    load V1k9; partmp = getObst(5); % From twoCircles, was without serpar
    par.obsts(1).serpar = partmp.obsts(1).serpar; par.obsts(2).serpar = partmp.obsts(2).serpar;
end
if 0
    load V1k128obst12swap; % Periodic orbit in other direction and swap back
    parSwap = par;
    par.obsts(2) = parSwap.obsts(3);
    par.obsts(3) = parSwap.obsts(2);
end

% Make sure tau_1 = 0 is in the middle
if 0
    l = length(par.obsts(1).colltau)/2;
    signal = [V1(l+1:end); V1(1:l)];
    collsignal = [(par.obsts(1).colltau(l+1:end)-1), par.obsts(1).colltau(1:l)];
else
    signal = V1;
    collsignal = par.obsts(1).colltau;
end
% figure; plot(collsignal, real(signal))


%% Calculate the phase
maxOrder = 5; 
[taus, c, a, ft, d] = seriesPhasePerOrbit(par, maxOrder);
tayPh = zeros(maxOrder, length(collsignal));
tayPh(1,:) = sum(d)*ones(size(collsignal));
for i = 1:maxOrder
    tayPh(i+1,:) = tayPh(i,:) + c(i,1)*(collsignal - taus(1)).^i;
end


%% Calculate phitilde: extract phase from angle
phitilde = zeros(size(signal));

closest = find(abs(collsignal -taus(1)) == min(abs(collsignal-taus(1))));
% phitilde(closest) = 1;
% phitilde(closest) = sum(d);
phitilde(closest) = angle(signal(closest))./par.k;

if 0
    for i = (closest-1):-1:1
        phitilde(i) = phitilde(i+1) +(angle(signal(i)) - angle(signal(i+1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i+1)))/2/pi);
        %     phitilde(i) = phitilde(i+1) +(angle(signal(i)) - angle(signal(i+1)))/par.k -2*pi/par.k*floor( abs(angle(signal(i)) -angle(signal(i+1)))/2/pi);
    end
    for i = (closest+1):length(signal) %length(V1)
        phitilde(i) = phitilde(i-1) +(angle(signal(i)) - angle(signal(i-1)))/par.k -2*pi/par.k*floor( (angle(signal(i)) -angle(signal(i-1)))/2/pi);
        %     phitilde(i) = phitilde(i-1) +(angle(signal(i)) - angle(signal(i-1)))/par.k -2*pi/par.k*floor( abs(angle(signal(i)) -angle(signal(i-1)))/2/pi);
    end
else
    multipl = 0;
    for i = (closest-1):-1:1
        phitilde(i) = angle(signal(i))./par.k + multipl*2*pi/par.k;
%         if phitilde(i) > (phitilde(i+1) +pi/par.k)
        if angle(signal(i)) > (angle(signal(i+1)) +pi)
            multipl = multipl -1;
            phitilde(i) = phitilde(i) - 2*pi/par.k;
%             phitilde(i) = phitilde(i) + 2*pi/par.k;
        elseif angle(signal(i)) < (angle(signal(i+1)) -pi)
            multipl = multipl +1;
            phitilde(i) = phitilde(i) + 2*pi/par.k;
        end
    end
    multipl = 0;
    for i = (closest+1):length(signal)
        phitilde(i) = angle(signal(i))./par.k + multipl*2*pi/par.k;
%         if phitilde(i) > (phitilde(i+1) +pi/par.k)
        if angle(signal(i)) > (angle(signal(i-1)) +pi)
            multipl = multipl -1;
            phitilde(i) = phitilde(i) - 2*pi/par.k;
%             phitilde(i) = phitilde(i) + 2*pi/par.k;
        elseif angle(signal(i)) < (angle(signal(i-1)) -pi)
            multipl = multipl +1;
            phitilde(i) = phitilde(i) + 2*pi/par.k;
        end
    end
end
% figure; plot(collsignal, phitilde);
% phitilde = sum(d) -(phitilde - phitilde(closest) );
phitilde = sum(d) +(phitilde - phitilde(closest) );


%% Plot phase

figure; 
phitPl = plot(collsignal, phitilde);
hold on;
% c0Pl = plot(collsignal, (phitilde- tayPh(1,:)')./phitilde);
c0Pl = plot(collsignal, tayPh(1,:));
c1Pl = plot(collsignal, tayPh(2,:));
c2Pl = plot(collsignal, tayPh(3,:));
c3Pl = plot(collsignal, tayPh(4,:));

% {'\tilde{\phi}$', '$(\tilde{\phi} -c_0)/\tilde{\phi}', '$(\tilde{\phi} -c_0 - c_1 [\tau_1 - \tau_1^*])/\tilde{\phi}',...
%     '$(\tilde{\phi} -\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i)/\tilde{\phi}', '$(\tilde{\phi} -\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i)/\tilde{\phi}'}, ...
% if 0 % gives error
legend([phitPl, c0Pl, c1Pl, c2Pl, c3Pl], {'$\tilde{\phi}$', '$c_0$', '$c_0 + c_1 [\tau_1 - \tau_1^*]$',...
    '$\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i$', '$\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i$'}, ...
    'interpreter', 'latex', 'FontSize', 20);

xlabel('\tau_{1,j}');
ylabel('Phase');
set(gca, 'FontSize', 20);

% return

%% Plot eigenvector
figure;
% V1plc8 = plot(collsignal, real(transpose(signal)./exp(1i*par.k*symbTay8(collsignal) )), 'g', 'LineWidth', 3);
% hold on;
% V1plc2 = plot(collsignal, real(transpose(signal)./exp(1i*par.k*symbTay2(collsignal))), 'r--', 'LineWidth', 2);
V1pl = plot(collsignal, real(transpose(signal)), 'b:');
% legend([V1pl, V1plc2, V1plc8], {'Re($V_{j,1})$', 'Re$\{V_{j,1}/\exp(ik[c_0+c_2\tau^2])\}$', 'Re$\{V_{j,1}/\exp(ik\sum_{i=0}^8 c_i \tau_1^i)\}$'},...
%     'interpreter', 'latex', 'FontSize', 20); %15);
xlabel('\tau_{1,j}');
ylabel('Mode');
set(gca, 'FontSize', 20);

%% Plot for second obstacle
sig2 = allV{2}(:,1);

phit2 = nan(size(sig2));
% Assume same collsignal
% clos2 = find(abs(collsignal -taus(2)) == min(abs(collsignal-taus(2))));
[~, clos2] = min(abs(collsignal-taus(2)) );

phit2(clos2) = angle(sig2(clos2))./par.k;

if 0
    for i = (clos2-1):-1:1
        phit2(i) = phit2(i+1) +(angle(sig2(i)) - angle(sig2(i+1)))/par.k -2*pi/par.k*floor( (angle(sig2(i)) -angle(sig2(i+1)))/2/pi);
    end
    for i = (clos2+1):length(phit2)
        phit2(i) = phit2(i-1) +(angle(sig2(i)) - angle(sig2(i-1)))/par.k -2*pi/par.k*floor( (angle(sig2(i)) -angle(sig2(i-1)))/2/pi);
    end
end
multipl = 0;
for i = (clos2-1):-1:1
    phit2(i) = angle(sig2(i))./par.k + multipl*2*pi/par.k;
    if angle(sig2(i)) > (angle(sig2(i+1)) +pi)
        multipl = multipl -1;
        phit2(i) = phit2(i) - 2*pi/par.k;
    elseif angle(sig2(i)) < (angle(sig2(i+1)) -pi)
        multipl = multipl +1;
        phit2(i) = phit2(i) + 2*pi/par.k;
    end
end
multipl = 0;
for i = (clos2+1):length(phit2)
    phit2(i) = angle(sig2(i))./par.k + multipl*2*pi/par.k;
    if angle(sig2(i)) > (angle(sig2(i-1)) +pi)
        multipl = multipl -1;
        phit2(i) = phit2(i) - 2*pi/par.k;
    elseif angle(sig2(i)) < (angle(sig2(i-1)) -pi)
        multipl = multipl +1;
        phit2(i) = phit2(i) + 2*pi/par.k;
    end
end
phit2 = sum(d) +(phit2 - phit2(clos2) );

tayPh2 = zeros(maxOrder, length(collsignal));
tayPh2(1,:) = sum(d)*ones(size(collsignal));
for i = 1:maxOrder-1
    tayPh2(i+1,:) = tayPh2(i,:) + c(i,2)*(collsignal - taus(2)).^i;
end

figure;
plot(collsignal, real(transpose(sig2)), 'b:');
xlabel('\tau_{2,j}');
ylabel('Mode in obst 2');
set(gca, 'FontSize', 20);

figure; 
phitPl = plot(collsignal, phit2);
hold on;
c0Pl = plot(collsignal, tayPh2(1,:));
c1Pl = plot(collsignal, tayPh2(2,:));
c2Pl = plot(collsignal, tayPh2(3,:));
c3Pl = plot(collsignal, tayPh2(4,:));

legend([phitPl, c0Pl, c1Pl, c2Pl, c3Pl], {'$\tilde{\phi}$', '$c_0$', '$c_0 + c_1 [\tau_1 - \tau_1^*]$',...
    '$\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i$', '$\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i$'}, ...
    'interpreter', 'latex', 'FontSize', 20);

xlabel('\tau_{1,j}');
ylabel('Phase for obst 2');
set(gca, 'FontSize', 20);

% convergence
figure; 
c0Pl = loglog(abs(collsignal-taus(2)), abs(phit2- tayPh2(1,:)')./phit2, 'b');
hold on;
c1Pl = plot(abs(collsignal-taus(2)), abs(phit2- tayPh2(2,:)')./phit2, 'r:');
c2Pl = plot(abs(collsignal-taus(2)), abs(phit2- tayPh2(3,:)')./phit2, 'g');
c3Pl = plot(abs(collsignal-taus(2)), abs(phit2- tayPh2(4,:)')./phit2, 'm-.');

legend([c0Pl, c1Pl, c2Pl, c3Pl], {'$|\tilde{\phi} -c_0|/\tilde{\phi}$', '$|\tilde{\phi} -c_0 - c_1 [\tau_1 - \tau_1^*]|/\tilde{\phi}$',...
    '$|\tilde{\phi} -\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i|/\tilde{\phi}$', '$|\tilde{\phi} -\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i|/\tilde{\phi}$'}, ...
    'interpreter', 'latex', 'FontSize', 20);

xlabel('|\tau_{2,j} - \tau_{2}^*|');
ylabel('Relative error for second obstacle');
set(gca, 'FontSize', 20);




%% Make plots
% figure; plot(collsignal, real(signal))

showTau = (0:5)/6;
perOrbit = nan(2,length(par.obsts));
ts = linspace(0,1,200)';

figure;
for moi = 1:length(par.obsts)
    topl = par.obsts(moi).par(ts');
    plot(topl(1,:), topl(2,:));
    hold on;
    h = text(mean(topl(1,:)) - 0.5*(max(topl(1,:))-min(topl(1,:)))+0.1, mean(topl(2,:)), ['Obst. ' num2str(moi)]);
    set(h,'FontSize',10);
    
    for si = 1:length(showTau)
        pt = par.obsts(moi).par(showTau(si));
        plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
%         h = text(pt(1)-0.3, pt(2)+0.1, ['$\tau=$ ' num2str(showTau(si))], 'interpreter', 'latex');
        h = text(pt(1) +0.04, pt(2), ['$\tau=$ ' num2str(showTau(si))], 'interpreter', 'latex');
        set(h,'FontSize',10, 'color', 'r');
    end
    perOrbit(:,moi) = par.obsts(moi).par(taus(moi));
end
plot([perOrbit(1,:) perOrbit(1,1)], [perOrbit(2,:) perOrbit(2,1)]);
set(h,'FontSize',10);
axis equal;


%% Plot relerr phase

figure; 
c0Pl = loglog(abs(collsignal-taus(1)), abs(phitilde- tayPh(1,:)')./phitilde, 'b');
hold on;
c1Pl = plot(abs(collsignal-taus(1)), abs(phitilde- tayPh(2,:)')./phitilde, 'r:');
c2Pl = plot(abs(collsignal-taus(1)), abs(phitilde- tayPh(3,:)')./phitilde, 'g');
c3Pl = plot(abs(collsignal-taus(1)), abs(phitilde- tayPh(4,:)')./phitilde, 'm-.');

legend([c0Pl, c1Pl, c2Pl, c3Pl], {'$|\tilde{\phi} -c_0|/\tilde{\phi}$', '$|\tilde{\phi} -c_0 - c_1 [\tau_1 - \tau_1^*]|/\tilde{\phi}$',...
    '$|\tilde{\phi} -\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i|/\tilde{\phi}$', '$|\tilde{\phi} -\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i|/\tilde{\phi}$'}, ...
    'interpreter', 'latex', 'FontSize', 20);

xlabel('|\tau_{1,j} - \tau_{1}^*|');
ylabel('Relative error');
set(gca, 'FontSize', 20);
