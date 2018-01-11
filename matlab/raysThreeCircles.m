% Plot the rays originating around the periodic orbit of three circles
clearvars
close all
clc
format longe
set(0,'DefaultFigureWindowStyle','docked');

mnr = 200;
tspl = linspace(0,1,mnr)';
par = getObst(12); % Three circles with changed dist
% par = getObst(5); % Two circles

[spTaus, c, a, ft] = seriesPhasePerOrbit(par, 3);
fss = 'Fontsize'; fs = 22;
lws = 'LineWidth'; lw = 5;
        
%% Start plot
perOrbitSPP = nan(2,length(par.obsts));
showTau = (0:5)/6;
cols = 'rgbmkcy';
lins = {'-', '--', ':', '-.'};
marks = {'v', '+', 'o', 'x', '*', 'h', 'd'};
figure;
hold on;
parametr = [];
shft = 0.3;
for moi = 1:length(par.obsts)
    topl = par.obsts(moi).par(tspl');
    plot(topl(1,:), topl(2,:));
    for si = 1:0 %length(showTau)
        pt = par.obsts(moi).par(showTau(si));
        plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
        h = text(pt(1) +0.04 -shft*(abs(showTau(si)-0.5) < 1/4), pt(2), ['$\tau=$ ' num2str(showTau(si), '%0.2f')], ...
            'interpreter', 'latex', 'FontSize', 16);
        set(h,'FontSize',20, 'color', 'r');
    end
    h = text(mean(topl(1,:)) - 0.5*(max(topl(1,:))-min(topl(1,:)))+0.1, mean(topl(2,:)), ['Obst. ' num2str(moi)], fss, fs);
    parametr = [parametr, topl];
    perOrbitSPP(:,moi) = par.obsts(moi).par(spTaus(moi));
end
plot([perOrbitSPP(1,:) perOrbitSPP(1,1)], [perOrbitSPP(2,:) perOrbitSPP(2,1)], 'b', 'LineWidth',3);
plot([0 0], [0.5 1.5], 'r', 'LineWidth', 3);
siz = max(parametr, [],2)-min(parametr,[],2);
mid = (max(parametr, [],2) +min(parametr,[],2) )/2;
axis equal


%% Reflections

% maxRefl = 9;
maxRefl = 6;
% maxRefl = 4;
% maxRefl = 3;
% maxRefl = 2;

% nray = 15;
nray = 8;

sel = @(x,y) x(y,:);
selco = @(x,y) x(:,y);
div = @(x) x(1,:)./x(2,:);
tol = 5e-4;

p2 = par.obsts(2).par(spTaus(2));
% Calculate point closest to p2, which is very near (but not equal) to the minimum of the phase
t1N = fminunc(@(ts) norm(p2 - [sin(2*pi*ts(1)); cos(2*pi*ts(1))]/2), 0, optimoptions('fminunc', 'Algorithm',...
    'quasi-newton', 'TolX', eps, 'TolFun', eps, 'OptimalityTolerance', eps, 'MaxFunctionEvaluations', 1000, 'Display', 'off'))
p1N = par.obsts(1).par(t1N);
% Do not plot because phase minimal for inf nr refl iso one  plot([p2(1), p1N(1)], [p2(2), p1N(2)], 'g', 'LineWidth',4);

active = 1:nray;
taus = nan(nray,maxRefl);
type = 0;
if type == 0 % Ray leaves orthogonally to the obstacle
%     taus(:,1) = linspace(-1,1,nray)*0.019646*(3-2*sqrt(2))^2; 
    taus(:,1) = linspace(0.0223,0.0228,nray);
    rat = div(par.obsts(1).normal(taus(:,1)'));
elseif type == 1 % Take one point with multiple directions iso multiple points with normal directions
    taus(:,1) = spTaus(1);
%     rat = cot(linspace(-pi/6,pi/3,nray));
    rat = cot(linspace(pi/2-0.03, pi/2+0.03,nray)); 
elseif type == 2 % Towards p2 at tau_{j+1}^* close to orth dir
    taus(:,1) = linspace(0,0.1,nray);
    rat = div(par.obsts(1).par(taus(:,1)') - repmat(par.obsts(2).par(spTaus(2)), [1, nray]) );
end

oi = 2;
pts = nan(2,nray,maxRefl);

pts(:,:,1) = par.obsts(1).par(taus(:,1)');

for ri = 2:maxRefl
    if 0 % vectorised fsolve
        [ts, fval] = fsolve(@(ts) deal(diag(arrayfun(@(ac) selco(sel(par.obsts(oi).serpar(ts(ac),2),2),2).*rat(ac) ...
            - selco(sel(par.obsts(oi).serpar(ts(ac),2),1),2), active)' ), ...
            pts(1,active,ri-1) -rat.*pts(2,active,ri-1) ...
            -sel(par.obsts(oi).par(ts),1) + rat.*sel(par.obsts(oi).par(ts),2) ), star(oi)*active.^0, optimoptions(@fsolve,'SpecifyObjectiveGradient',true) );
        active = active(abs(fval) < tol)
        taus(active,ri) = ts(abs(fval) < tol);
    else
        toRem = [];
        for ai = 1:length(active)
            if 0 % Use gradient
                [ts, fval] = fsolve(@(ts) deal( pts(1, active(ai), ri-1) -rat(ai).*pts(2, active(ai), ri-1) ...
                    -sel(par.obsts(oi).par(ts),1) + rat(ai).*sel(par.obsts(oi).par(ts),2)    ,     selco(sel(par.obsts(oi).serpar(ts,2), 2), 2).*rat(ai) ...
                    - selco(sel(par.obsts(oi).serpar(ts, 2), 1), 2)  ), star(oi), optimoptions(@fsolve,'SpecifyObjectiveGradient',true, 'Display','off') );
            else
                [ts, fval] = fsolve(@(ts) pts(1, active(ai), ri-1) -sel(par.obsts(oi).par(ts),1)   ...
                    + rat(ai).*(sel(par.obsts(oi).par(ts),2) -pts(2, active(ai), ri-1)) , spTaus(oi), ...
                    optimoptions(@fsolve, 'Display','off')  );
            end
%             display([num2str(active(ai)) ' = active(ai), fval = ' num2str([ts,fval])])
            if abs(fval) < tol
                taus(active(ai),ri) = ts;
            else
                toRem = [toRem, ai];
            end
        end
        active(toRem) = [];
        display(taus(active,1)');
    end
    
    pts(:, active ,ri) = par.obsts(oi).par(taus(active,ri)');
    rat = cot(2*acot(div(par.obsts(oi).normal(taus(active,ri)'))) ... 
        - acot(div(pts(:,active,ri-1) - pts(:,active,ri) )) );
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(squeeze(pts(1, :, (ri-1):ri))', squeeze(pts(2, :, (ri-1):ri))'); % Not pts(1,active to make sure same colour
    oi = oi + 1;
    if oi > length(par.obsts), oi = 1; end
end
set(gca, fss, fs);
