% Make convergence plots for the phase of the periodic orbit of general obstacles.
% addpath('/dir/to/chebfun')

%% Loading saved eigenvector
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');
% Loads V1l, allV and par
% load V1k128obst12; % Three circles
% load V1k128obst13; % Nearconvex and ellipse
% load V1k512obst13; % Nearconvex and ellipse at higher frequency
load V1k128obst14; % Nearconvex, ellipse and near-inclusion
% load V1k128obst5.mat; % Two circles
% load V1k128obst2nsEll; % Two nonsymmetric ellipses
% load V1k128obst3nsCirc; % Three nonsymmetric circles
% load V1k128obstcircFFtEl; % Circle with FFT and ellipse
% load V1k128obstcircWoFFtElswap; % Circle without FFT and ellipse
% load V1k128obstnincNcEll; % Near-inclusion near convex part, near-convex and ellipse
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
% maxOrder = 5; 
% maxOrder = 20;
maxOrder = 6;
[taus, c, a, ft, d] = seriesPhasePerOrbit(par, maxOrder);
tayPh = zeros(maxOrder, length(collsignal));
tayPh(1,:) = sum(d)*ones(size(collsignal));
for i = 1:maxOrder
    tayPh(i+1,:) = tayPh(i,:) + c(i,1)*(collsignal - taus(1)).^i;
end

if 0 % Do for any obstacle below
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
end

%% Plot for some obstacle with index testo
% testo = 1;
% testo = 2;
% testo = 3;
for testo = 1:length(par.obsts)
% sig2 = allV{2}(:,1);
sig2 = allV{testo}(:,1);

phit2 = nan(size(sig2));
% Assume same collsignal
% clos2 = find(abs(collsignal -taus(2)) == min(abs(collsignal-taus(2))));
% [~, clos2] = min(abs(collsignal-taus(2)) );
[~, clos2] = min(abs(collsignal-taus(testo)) );

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
% Use periodized version of collsignal in series expansion
collser = collsignal + round(taus(testo) - collsignal);
for i = 1:maxOrder-1
%     tayPh2(i+1,:) = tayPh2(i,:) + c(i,2)*(collsignal - taus(2)).^i;
%     tayPh2(i+1,:) = tayPh2(i,:) + c(i,testo)*(collsignal - taus(testo)).^i;
    tayPh2(i+1,:) = tayPh2(i,:) + c(i,testo)*(collser - taus(testo)).^i;
end

% marks = {'b', 'r:', 'g--', 'k-.', 'c-', 'm:', 'y'};
% marks = {'b', 'r:', 'g--', 'k-.', 'c-', 'm-.', 'y'};
marks = {'b', 'r:', 'g--', 'k-.', 'c-'};
ords = unique(round((1:length(marks))*maxOrder/length(marks)));
mg = 'b:';
shft = 0.32; % For nincNcEll
if abs(taus(1) -6.174348705381557e-01) + abs(taus(2) -6.747611126692694e-02) < 1e-5 % then ncEll
    marks = {'k-.', 'r--', 'g-.', 'c', 'b--'};
    ords = [3 6];
    shft = 0.3;
%     shft = 0.16;
elseif abs(taus(1) - 7.276256269263204e-01) + abs(taus(2) - 2.734466122469346e-02) ...
        + abs(taus(min(3,length(par.obsts))) - 3.090147854459014e-01) < 1e-5 %obst 14
    marks = {'r:', 'g--', 'k-.', 'c-'};
    shft = 0.4;
    ords = [4 6];
    mg = 'b';
%     shft = 0.24;
end
figure;
% V1plPh = cell(maxOrder,1);
% V1plPh = gobjects(maxOrder+1,1);
% legs = cell(maxOrder+1,1);
V1plPh = gobjects(length(ords) +1,1);
legs = cell(length(ords) +1,1);
% V1plPh(1) = plot(collsignal, real(transpose(sig2)), 'b:');
V1plPh(1) = plot(collsignal, real(transpose(sig2)), mg);
hold on;
legs{1} = 'Re($V_{j,1})$';
% for ord = 1:maxOrder %maxOrder:-1:1
%     V1plPh(ord+1) = plot(collsignal, real(transpose(sig2)./exp(1i*par.k*tayPh2(ord,:) )), marks{ord}, 'LineWidth', 3);
%     legs{ord+1} = ['Re$\{V_{j,1}/\exp(ik\sum_{i=0}^{' num2str(ord-1) '} c_i [\tau - \tau^*]^i)\}$'];
% for oi = 1:length(ords)
for oi = length(ords):-1:1
    ord = ords(oi);
    V1plPh(oi+1) = plot(collsignal, real(transpose(sig2)./exp(1i*par.k*tayPh2(ord,:) )), marks{oi}, ...
        'LineWidth', oi/length(ords)*2 +1);
    legs{oi+1} = ['Re$\{V_{j,1}/\exp(ik\sum_{i=0}^{' num2str(ord-1) '} c_i [\tau - \tau^*]^i)\}$'];
end
% V1plPh(1) = plot(collsignal, real(transpose(sig2)), 'b:');
% legs{1} = 'Re($V_{j,1})$';
% V1plcmax = plot(collsignal, real(transpose(sig2)./exp(1i*par.k*tayPh2(maxOrder,:) )), 'g', 'LineWidth', 3);
% hold on;
% ho = round(maxOrder/2);
% V1plcha = plot(collsignal, real(transpose(sig2)./exp(1i*par.k*tayPh2(ho,:) )), 'r--', 'LineWidth', 2);
% V1pl = plot(collsignal, real(transpose(sig2)), 'b:');
% legend([V1pl, V1plcha, V1plcmax], {'Re($V_{j,1})$', ['Re$\{V_{j,1}/\exp(ik\sum_{i=0}^{' num2str(ho) '} c_i [\tau - \tau^*]^i)\}$'], ...
%     ['Re$\{V_{j,1}/\exp(ik\sum_{i=0}^{' num2str(maxOrder) '} c_i [\tau - \tau^*]^i)\}$']},...
%     'interpreter', 'latex', 'FontSize', 20); 

legend(V1plPh, legs, 'interpreter', 'latex', 'FontSize', 20, 'location', 'best');
xlabel(['\tau_{' num2str(testo) ',j}']);
% xlabel('\tau_{2,j}');
ylabel(['Mode in obstacle ' num2str(testo)]);
% ylabel('Mode in obst 2');
% set(gca, 'FontSize', 20);
set(gca, 'FontSize', 26);


% Plot the convergence of the phase
if abs(taus(1) -6.174348705381557e-01) + abs(taus(2) -6.747611126692694e-02) < 1e-5 % then ncEll
    ords = [1 3 4 5 6];
elseif abs(taus(1) - 7.276256269263204e-01) + abs(taus(2) - 2.734466122469346e-02) ...
        + abs(taus(min(3,length(par.obsts))) - 3.090147854459014e-01) < 1e-5 %obst 14
    ords = [1 4 6];
    marks = {'b', 'r:', 'g--', 'k-.', 'c-'};
end
% marks = {'b', 'r:', 'g--', 'k-.', 'c-', 'm:', 'y'};
figure;
% cvgs = gobjects(maxOrder,1);
% legs = cell(maxOrder,1);
cvgs = gobjects(length(ords), 1);
legs = cell(length(ords), 1);
% for ord = 1:maxOrder
%     cvgs(ord) = loglog(abs(collsignal-taus(testo)), abs(phit2- tayPh2(ord,:)')./phit2, marks{ord});
%     legs{ord} = ['$|\tilde{\phi} -\sum_{i=0}^{' num2str(ord-1) '} c_i [\tau - \tau^*]^i|/\tilde{\phi}$'];
for oi = length(ords):-1:1
    ord = ords(oi);
    cvgs(oi) = loglog(abs(collsignal-taus(testo)), abs(phit2- tayPh2(ord,:)')./phit2, marks{oi},'LineWidth', oi/length(ords)*2 +1);
    hold on;
    legs{oi} = ['$|\tilde{\phi} -\sum_{i=0}^{' num2str(ord-1) '} c_i [\tau - \tau^*]^i|/\tilde{\phi}$'];
end

legend(cvgs, legs, 'interpreter', 'latex', 'FontSize', 20, 'location', 'best');
xlabel(['|\tau_{' num2str(testo) ',j} - \tau_{' num2str(testo) '}^*|']);
ylabel(['Relative error for obstacle ' num2str(testo)]);
set(gca, 'FontSize', 26);
% set(gca, 'FontSize', 20);

other = testo-1;
if other == 0
    other = length(par.obsts);
end
minDist = fminsearch(@(t) norm(par.obsts(other).par(taus(other)) -par.obsts(testo).par(t) ), taus(testo));
[mi, ix] = min(phit2);
[collsignal(ix), minDist, (collsignal(ix) -minDist), (collsignal(ix) -collsignal(ix-1))]
end

%% Old plots

if 0
figure; 
phitPl = plot(collsignal, phit2);
hold on;
c0Pl = plot(collsignal, tayPh2(1,:));
c1Pl = plot(collsignal, tayPh2(2,:));
c2Pl = plot(collsignal, tayPh2(3,:));
c3Pl = plot(collsignal, tayPh2(4,:));

% legend([phitPl, c0Pl, c1Pl, c2Pl, c3Pl], {'$\tilde{\phi}$', '$c_0$', '$c_0 + c_1 [\tau_1 - \tau_1^*]$',...
%     '$\sum_{i=0}^2 c_i [\tau_1 - \tau_1^*]^i$', '$\sum_{i=0}^3 c_i [\tau_1 - \tau_1^*]^i$'}, ...
legend([phitPl, c0Pl, c1Pl, c2Pl, c3Pl], {'$\tilde{\phi}$', '$c_0$', '$c_0 + c_1 [\tau - \tau^*]$',...
    '$\sum_{i=0}^2 c_i [\tau - \tau^*]^i$', '$\sum_{i=0}^3 c_i [\tau - \tau^*]^i$'}, ...
    'interpreter', 'latex', 'FontSize', 20);
xlabel(['\tau_{' num2str(testo) ',j}']);
ylabel(['Phase for obst ' num2str(testo)]);
set(gca, 'FontSize', 20);
% end

% convergence
figure; 
c0Pl = loglog(abs(collsignal-taus(testo)), abs(phit2- tayPh2(1,:)')./phit2, 'b');
% c0Pl = loglog(abs(collsignal-taus(2)), abs(phit2- tayPh2(1,:)')./phit2, 'b');
hold on;
c1Pl = plot(abs(collsignal-taus(testo)), abs(phit2- tayPh2(2,:)')./phit2, 'r:');
c2Pl = plot(abs(collsignal-taus(testo)), abs(phit2- tayPh2(3,:)')./phit2, 'g');
c3Pl = plot(abs(collsignal-taus(testo)), abs(phit2- tayPh2(4,:)')./phit2, 'm-.');

legend([c0Pl, c1Pl, c2Pl, c3Pl], {'$|\tilde{\phi} -c_0|/\tilde{\phi}$', '$|\tilde{\phi} -c_0 - c_1 [\tau - \tau^*]|/\tilde{\phi}$',...
    '$|\tilde{\phi} -\sum_{i=0}^2 c_i [\tau - \tau^*]^i|/\tilde{\phi}$', '$|\tilde{\phi} -\sum_{i=0}^3 c_i [\tau - \tau^*]^i|/\tilde{\phi}$'}, ...
    'interpreter', 'latex', 'FontSize', 20);

% xlabel('|\tau_{2,j} - \tau_{2}^*|');
% ylabel('Relative error for second obstacle');
xlabel(['|\tau_{' num2str(testo) ',j} - \tau_{' num2str(testo) '}^*|']);
ylabel(['Relative error for obstacle ' num2str(testo)]);
set(gca, 'FontSize', 20);

end

%% Make plots
% figure; plot(collsignal, real(signal))

showTau = (0:5)/6;
perOrbit = nan(2,length(par.obsts));
ts = linspace(0,1,200)';

% shft = 0.32; % For nincNcEll
% if abs(taus(1) -6.174348705381557e-01) + abs(taus(2) -6.747611126692694e-02) < 1e-5 % then ncEll
%     shft = 0.16;
% end

figure;
for moi = 1:length(par.obsts)
    topl = par.obsts(moi).par(ts');
    plot(topl(1,:), topl(2,:));
    hold on;
%     h = text(mean(topl(1,:)) - 0.5*(max(topl(1,:))-min(topl(1,:)))+0.1, mean(topl(2,:)), ['Obst. ' num2str(moi)]);
    h = text(mean(topl(1,:))- 0.4*(max(topl(1,:))-min(topl(1,:))), ...
        mean(topl(2,:)) - 0.12*(max(topl(2,:))-min(topl(2,:))), ['Obst. ' num2str(moi)]);
    set(h,'FontSize',20);
    
    for si = 1:length(showTau)
        pt = par.obsts(moi).par(showTau(si));
        plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
%         h = text(pt(1)-0.3, pt(2)+0.1, ['$\tau=$ ' num2str(showTau(si))], 'interpreter', 'latex');
%         h = text(pt(1) +0.04 -0.32*(abs(showTau(si)-0.5) < 1/4), pt(2), ['$\tau=$ ' num2str(showTau(si), '%0.2f')], 'interpreter', 'latex');
        h = text(pt(1) +0.04 -shft*(abs(showTau(si)-0.5) < 1/4), pt(2), ['$\tau=$ ' num2str(showTau(si), '%0.2f')], ...
            'interpreter', 'latex', 'FontSize', 16);
%         h = text(pt(1) +0.04, pt(2), ['$\tau=$ ' num2str(showTau(si))], 'interpreter', 'latex');
%         h = text(pt(1) +0.04 - 0.5*(abs(showTau(moi)-0.5) < 1/4), pt(2), ['$\tau=$ ' num2str(showTau(si))], 'interpreter', 'latex');
        set(h,'FontSize',20, 'color', 'r');
    end
    perOrbit(:,moi) = par.obsts(moi).par(taus(moi));
end
plot([perOrbit(1,:) perOrbit(1,1)], [perOrbit(2,:) perOrbit(2,1)]);
% quiver( (perOrbit(1,1:end-1) + perOrbit(1,2:end))/2, (perOrbit(2,1:end-1) + perOrbit(2,2:end))/2,...
%     (perOrbit(1,2:end) - perOrbit(1,1:end-1))/6, (perOrbit(2,2:end) - perOrbit(2,1:end-1))/6);
% quiver( (perOrbit(1,1) + perOrbit(1,end))/2, (perOrbit(2,1) + perOrbit(2,end))/2, ...
%      (perOrbit(1,1)- perOrbit(1,end))/6, (perOrbit(2,1)- perOrbit(2,end))/6);
 quiver( [(perOrbit(1,1:end-1) + perOrbit(1,2:end)), (perOrbit(1,1) + perOrbit(1,end))]/2, ...
    [(perOrbit(2,1:end-1) + perOrbit(2,2:end)), (perOrbit(2,1) + perOrbit(2,end))]/2,...
    [(perOrbit(1,2:end) - perOrbit(1,1:end-1)), (perOrbit(1,1)- perOrbit(1,end))]/16, ...
    [(perOrbit(2,2:end) - perOrbit(2,1:end-1)), (perOrbit(2,1)- perOrbit(2,end))]/16);
% set(h,'FontSize',10);
axis equal;
set(gca, 'FontSize', 20);

if shft == 0.3
    yticks([0.4, 0.6, 0.8]);
elseif shft == 0.4
    yticks([-0.6,  -0.4, 0.4, 0.6, 0.8, 1]);
end

if 0
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
end

