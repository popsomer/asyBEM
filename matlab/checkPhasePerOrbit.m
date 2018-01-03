% Make convergence plots for the phase of the periodic orbit of general obstacles.
% First run makeVb1.m with the desired obstacle and wavenumber, then add the resulting 
% .mat below and possibly adjust some plotting parameters.

%% Loading saved eigenvector: .mat contains V1l, allV and par
clearvars
close all
format longe
set(0,'DefaultFigureWindowStyle','docked');

% load V1k512obst13; % Nearconvex and ellipse at higher frequency: run makeVb1.m for this
load V1k128obst14; % Nearconvex, ellipse and near-inclusion

signal = V1;
collsignal = par.obsts(1).colltau;


%% Calculate the phase series coefficients
maxOrder = 6;
[taus, c, a, ft, d] = seriesPhasePerOrbit(par, maxOrder);


%% Plot for some obstacle with index testo

for testo = 1:length(par.obsts)
    sig2 = allV{testo}(:,1);
    
    phit2 = nan(size(sig2));
    % Assume same collsignal
    % clos2 = find(abs(collsignal -taus(2)) == min(abs(collsignal-taus(2))));
    % [~, clos2] = min(abs(collsignal-taus(2)) );
    [~, clos2] = min(abs(collsignal-taus(testo)) );
    
    phit2(clos2) = angle(sig2(clos2))./par.k;
    
    % Extract the phase from the eigenvector
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
        tayPh2(i+1,:) = tayPh2(i,:) + c(testo,i)*(collser - taus(testo)).^i;
    end
    
    % Select plotting parameters
    marks = {'b', 'r:', 'g--', 'k-.', 'c-'};
    ords = unique(round((1:length(marks))*maxOrder/length(marks)));
    mg = 'b:';
    shft = 0.32; % For nincNcEll
    if abs(taus(1) -6.174348705381557e-01) + abs(taus(2) -6.747611126692694e-02) < 1e-5 % then ncEll
        marks = {'k-.', 'r--', 'g-.', 'c', 'b--'};
        ords = [3 6];
        shft = 0.3;
    elseif abs(taus(1) - 7.276256269263204e-01) + abs(taus(2) - 2.734466122469346e-02) ...
            + abs(taus(min(3,length(par.obsts))) - 3.090147854459014e-01) < 1e-5 %obst 14
        marks = {'r:', 'g--', 'k-.', 'c-'};
        shft = 0.4;
        ords = [4 6];
        mg = 'b';
    end
    
    figure;
    V1plPh = gobjects(length(ords) +1,1);
    legs = cell(length(ords) +1,1);
    V1plPh(1) = plot(collsignal, real(transpose(sig2)), mg);
    hold on;
    legs{1} = ['Re$\{\tilde{V}_{' num2str(testo) '}(\tau_{' num2str(testo) '})$'];
    for oi = length(ords):-1:1
        ord = ords(oi);
        V1plPh(oi+1) = plot(collsignal, real(transpose(sig2)./exp(1i*par.k*tayPh2(ord,:) )), marks{oi}, ...
            'LineWidth', oi/length(ords)*2 +1);
        legs{oi+1} = ['Re$\{\tilde{V}_{' num2str(testo) '}(\tau_{' num2str(testo) '})/\exp(ik\sum_{i=0}^{' num2str(ord-1) '} c_{' num2str(testo) ...
            ',i} [\tau_{' num2str(testo) '} - \tau_{' num2str(testo) '}^*]^i)\}$'];
    end
    
    legend(V1plPh, legs, 'interpreter', 'latex', 'FontSize', 20, 'location', 'best');
    xlabel(['\tau_{' num2str(testo) '}']);
    ylabel(['Mode in obstacle ' num2str(testo)]);
    set(gca, 'FontSize', 26);
    
    
    % Plot the convergence of the phase
    if abs(taus(1) -6.174348705381557e-01) + abs(taus(2) -6.747611126692694e-02) < 1e-5 % then ncEll
        ords = [1 3 4 5 6];
    elseif abs(taus(1) - 7.276256269263204e-01) + abs(taus(2) - 2.734466122469346e-02) ...
            + abs(taus(min(3,length(par.obsts))) - 3.090147854459014e-01) < 1e-5 %obst 14
        ords = [1 4 6];
        marks = {'b', 'r:', 'g--', 'k-.', 'c-'};
    end
    figure;
    cvgs = gobjects(length(ords), 1);
    legs = cell(length(ords), 1);
    for oi = length(ords):-1:1
        ord = ords(oi);
        cvgs(oi) = loglog(abs(collsignal-taus(testo)), abs(phit2- tayPh2(ord,:)')./phit2, marks{oi},'LineWidth', oi/length(ords)*2 +1);
        hold on;
        legs{oi} = ['$|\tilde{\phi}_{' num2str(testo) '}(\tau_{' num2str(testo) '}) -\sum_{i=0}^{' num2str(ord-1) '} c_{' ...
            num2str(testo) ',i} [\tau_{' num2str(testo) '} - \tau_{' num2str(testo) '}^*]^i|/\tilde{\phi}_{' num2str(testo) '}(\tau_{' num2str(testo) '})$'];
    end
    
    legend(cvgs, legs, 'interpreter', 'latex', 'FontSize', 20, 'location', 'best');
    xlabel(['|\tau_{' num2str(testo) '} - \tau_{' num2str(testo) '}^*|']);
    ylabel(['Relative error for obstacle ' num2str(testo)]);
    xticks(10.^(-4:0));
    yticks(10.^(-8:2:2));
    set(gca, 'FontSize', 26);
    
    other = testo-1;
    if other == 0
        other = length(par.obsts);
    end
    minDist = fminsearch(@(t) norm(par.obsts(other).par(taus(other)) -par.obsts(testo).par(t) ), taus(testo));
    [mi, ix] = min(phit2);
    [collsignal(ix), minDist, (collsignal(ix) -minDist), (collsignal(ix) -collsignal(ix-1))]
end


%% Plot the obstacles

showTau = (0:5)/6;
perOrbit = nan(2,length(par.obsts));
ts = linspace(0,1,200)';

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
        h = text(pt(1) +0.04 -shft*(abs(showTau(si)-0.5) < 1/4), pt(2), ['$\tau=$ ' num2str(showTau(si), '%0.2f')], ...
            'interpreter', 'latex', 'FontSize', 16);
        set(h,'FontSize',20, 'color', 'r');
    end
    perOrbit(:,moi) = par.obsts(moi).par(taus(moi));
end
plot([perOrbit(1,:) perOrbit(1,1)], [perOrbit(2,:) perOrbit(2,1)]);
 quiver( [(perOrbit(1,1:end-1) + perOrbit(1,2:end)), (perOrbit(1,1) + perOrbit(1,end))]/2, ...
    [(perOrbit(2,1:end-1) + perOrbit(2,2:end)), (perOrbit(2,1) + perOrbit(2,end))]/2,...
    [(perOrbit(1,2:end) - perOrbit(1,1:end-1)), (perOrbit(1,1)- perOrbit(1,end))]/16, ...
    [(perOrbit(2,2:end) - perOrbit(2,1:end-1)), (perOrbit(2,1)- perOrbit(2,end))]/16);
axis equal;
set(gca, 'FontSize', 20);

if shft == 0.3
    yticks([0.4, 0.6, 0.8]);
elseif shft == 0.4
    yticks([-0.6,  -0.4, 0.4, 0.6, 0.8, 1]);
end
