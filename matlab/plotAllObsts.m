% Plot all obstacles.
%% Initialising
clearvars
close all
clc
format longe
set(0,'DefaultFigureWindowStyle','docked');


%% All obstacles
mnr = 200;
ts = linspace(0,1,mnr)';
showTau = (0:5)/6;
labels = {'Circle', 'Ellipse', 'Near-inclusion', 'Almost convex', ...
    'Two circles', 'Three ellipses', 'Circle and near-inclusion', 'Nonconvex polygon'};
las = {'\xi=4%, N/k=6', '\xi=4%, N/k=6', '\xi=1%, N/k=10', '\xi=4%, N/k=7', ...
    '\xi=0.3%, N/k/2=10', '\xi=0.3%, N/k/3=10', '\xi=0.1%, N/k/2=10', '\xi=4%, N/k=9'};
hormove = [0.3, 0.05,0, -0.2  , -1, -0.2, -0.3, 0];
vertmove = [0, -0.1, 0.4, -0.1  , 0.5 , 1.1, 0.2, 0];

cols = 'rgbmkcy';
marks = {'v', '+', 'o', 'x', '*', 'h', 'd'};
figure; 
for oi = 3:8
    par = getObst(oi);
    subplot(2,3,oi-2);
    if isfield(par,'obsts')
        parametr = [];
        for moi = 1:length(par.obsts)
            topl = par.obsts(moi).par(ts');
            plot(topl(1,:), topl(2,:));
            hold on;
            h = text(mean(topl(1,:)) - 0.5*(max(topl(1,:))-min(topl(1,:)))+0.1, mean(topl(2,:)), ['Obst. ' num2str(moi)]);
            set(h,'FontSize',10);
            parametr = [parametr, topl];
            topl = par.obsts(moi).par(showTau);
            if (oi == 6) && (moi == 3)
                pt = par.obsts(moi).par(0.5);
                plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
                h = text(pt(1)-0.3, pt(2)+0.1, '$\mathbf{p}$', 'interpreter', 'latex');
                set(h,'FontSize',10, 'color', 'r');
            end
        end
    else
        topl = par.par(showTau);
        parametr = par.par(ts');
        plot(parametr(1,:), parametr(2,:));
        if oi == 3
            hold on;
            pt = par.par(0.35);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1)+0.05, pt(2)+0.04, '$\mathbf{q}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
            
            pt = par.par(0.5);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1)+0.02, pt(2)-0.02, '$\mathbf{r}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
            
            pt = par.par(0.65);
            plot(pt(1), pt(2), 'ro','MarkerFaceColor','r', 'MarkerSize',5);
            h = text(pt(1), pt(2)-0.07, '$\mathbf{s}$', 'interpreter', 'latex');
            set(h,'FontSize',10, 'color', 'r');
        end
    end
    siz = max(parametr, [],2)-min(parametr,[],2);
    mid = (max(parametr, [],2) +min(parametr,[],2) )/2;
    title({labels{oi}; ['\rm{' las{oi} '}']});
    set(h,'FontSize',10);
    axis([mid(1)-0.6*max(siz), mid(1)+0.6*max(siz), mid(2)-0.6*max(siz), mid(2)+0.6*max(siz)]);
end
