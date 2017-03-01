% Plot the validation results.
% Input
%   v		- Validation data: all fields are optional and see validate.m for their description.
%   [A		- Matrix of which to plot the sparsity structure, or empty if not]
%   [idxs   - If idxs is a vector, it gives the range of indices in v to plot, probably referring to the tested wavenumbers for one obstacle. 
%             If it is a cell, we print a LaTeX table of all results using the corresponding names when A=0 or figures of them when A=1 ]
function plotVal(v,A,idxs)
fss = 'Fontsize'; fs = 22;
lws = 'LineWidth'; lw = 5;
mss = 'MarkerSize'; ms = 15;
if ~exist('idxs','var')
    if isfield(v,'ks')
        idxs = 1:length(v.ks);
    else
        idxs = [];
    end
end
l = {'+', 'x', 'o', 'v', '*', 'h', 'd'}; 
ll = {'-','--',':', '-.','-','--',':','-.'};
c = 'brgcmy';

solMeth = {'$A$\\textbackslash$b$','$\\tilde{A}$\\textbackslash$b$','GM $A$\\textbackslash$b$','GM $\\tilde{A}$\\textbackslash$b$',...
        'GM $(\\tilde{A}^{-1}A)$\\textbackslash$(\\tilde{A}^{-1}b)$', 'GM $A$\\textbackslash$b$ with $x_0=\\tilde{x}$'};

if iscell(idxs) && (A == 0) % Print a LaTeX table of results for all obstacles.
    no = length(idxs);
    kl = length(v.ks);
    nbs = 2;
    if isfield(v,'mti')
        nbs = 2+v.mti;
    end
    fprintf('\\begin{table}[h] \n')
    fprintf('\\centering \n')
    fprintf(['\\hspace*{0.0cm} \\begin{tabular}{l|' repmat('c',1,kl) '} \n'])
    str =  sprintf('%u & ', v.ks);
    fprintf(['k & ' str(1:end-2) '\\\\ \n']);
    fprintf('\\hline \n')
    if isfield(v, 'errSol')
        for oi = 1:no
            for ns = 2:nbs
                str = regexprep(sprintf('%1.2e & ',v.errSol(kl*(oi-1)+(1:kl),ns)'),'e-0','e-');
                fprintf([idxs{oi} ' $||c-($ ' solMeth{ns} ' $)||/||c||$  & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'perc')
        for oi = 1:no
            str = sprintf('%2.1f & ', 100*v.perc(kl*(oi-1)+(1:kl),2)');
            fprintf([idxs{oi} ' \\%% nnz  & ' str(1:end-2) '\\\\ \n']);
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'compresErr')
        for oi = 1:no
            str = regexprep(sprintf('%1.2e & ', v.compresErr(kl*(oi-1)+(1:kl),2)'),'e-0','e-');
            fprintf([idxs{oi} ' $||\\tilde{A}c-b||/||b||$  & ' str(1:end-2) '\\\\ \n']);
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'errBCavm')
        for oi = 1:no
            for ns = 1:nbs
                str = regexprep(sprintf('%1.2e & ', v.errBCavm(kl*(oi-1)+(1:kl),ns)'),'e-0','e-');
                fprintf([idxs{oi} ' err BC ' solMeth{ns} ' & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    fprintf('\\end{tabular} \n')
    fprintf('\\caption{part 1} \n')
    fprintf('\\end{table} \n')
    
    fprintf('\\begin{table}[h] \n')
    fprintf('\\centering \n')
    fprintf(['\\hspace*{0.0cm} \\begin{tabular}{l|' repmat('c',1,kl) '} \n'])
    str =  sprintf('%u & ', v.ks);
    fprintf(['k & ' str(1:end-2) '\\\\ \n']);
    fprintf('\\hline \n')
    if isfield(v, 'errInt')
        for oi = 1:no
            for ns = 1:nbs
                str = regexprep(sprintf('%1.2e & ', v.errInt(kl*(oi-1)+(1:kl),ns)'),'e-0','e-');
                fprintf([idxs{oi} ' interior field ' solMeth{ns} ' & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'conds')
        for oi = 1:no
            str = regexprep(sprintf('%1.2e & ', v.conds(kl*(oi-1)+(1:kl),1)'),'e\+0','e');
            fprintf([idxs{oi} ' cond$(A)$  & ' str(1:end-2) '\\\\ \n']);
            str = regexprep(sprintf('%1.2e & ', v.conds(kl*(oi-1)+(1:kl),2)'),'e\+0','e');
            fprintf([idxs{oi} ' cond$(\\tilde{A})$  & ' str(1:end-2) '\\\\ \n']);
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'iterGm')
        for oi = 1:no
            for ng = 1:v.nbGm
                str = sprintf('%d & ',v.iterGm(kl*(oi-1)+(1:kl),ng)');
                fprintf([idxs{oi} ' \\# GMRES it. ' solMeth{2+ng} ' & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'nbIter')
        for oi = 1:no
            for ng = 1:v.mti
                str = sprintf('%d & ',v.nbIter(kl*(oi-1)+(1:kl),ng)');
                fprintf([idxs{oi} ' \\# GMRES it. ' solMeth{2+ng} ' & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    fprintf('\\end{tabular} \n')
    fprintf('\\caption{part 2} \n')
    fprintf('\\end{table} \n')
    
    fprintf('\\begin{table}[h] \n')
    fprintf('\\centering \n')
    fprintf(['\\hspace*{0.0cm} \\begin{tabular}{l|' repmat('c',1,kl) '} \n'])
    str =  sprintf('%u & ', v.ks);
    fprintf(['k & ' str(1:end-2) '\\\\ \n']);
    fprintf('\\hline \n')
    if isfield(v, 'timeSol')
        for oi = 1:no
            for ns = 1:nbs
                str = regexprep(regexprep(sprintf('%1.2e & ',v.timeSol(kl*(oi-1)+(1:kl),ns)') ,'e-0','e-'), 'e\+0', 'e');
                fprintf([idxs{oi} ' time ' solMeth{ns} ' (s) & ' str(1:end-2) '\\\\ \n']);
            end
        end
        fprintf('\\hline \n');
    end
    if isfield(v, 'timeA')
        nam = {'$A$','$\\tilde{A}$', 'Total (including validation) ','$r_{m,n}$'};
        for oi = 1:no
            for ns = 1:size(v.timeA,2)
                str = regexprep(regexprep(sprintf('%1.2e & ',v.timeA(kl*(oi-1)+(1:kl),ns)') ,'e-0','e-'), 'e\+0', 'e');
                fprintf([idxs{oi} ' time ' nam{ns} ' (s) & ' str(1:end-2) '\\\\ \n']);
            end
        end
    end
    fprintf('\\end{tabular} \n')
    fprintf('\\caption{.} \n')
    fprintf('\\label{} \n')
    fprintf('\\end{table} \n')
    return
elseif iscell(idxs) % Plot results for all obstacles.
    nrObsts = size(v.conds,1)/length(v.ks);
    
    figure;
    for oi = 1:nrObsts
        loglog(v.ks, v.errSol(length(v.ks)*(oi-1)+(1:length(v.ks)),2), ll{oi}, lws, lw);
        hold on;
    end
    xlabel('k',fss,fs); legend(idxs, 'interpreter','latex');
	ylabel('$|| \tilde{c}-$c$ ||/|| $c$ ||$','interpreter','latex',fss,fs); set(gca,fss,fs);
    
    mylegds = cell(2*nrObsts,1);
    figure;
    for oi = 1:nrObsts
        loglog(v.ks, v.errBCavm(length(v.ks)*(oi-1)+(1:length(v.ks)),:), ll{oi}, lws, lw);
        mylegds{2*(oi-1)+1} = [idxs{oi} ' $A$\textbackslash$b$'];
        mylegds{2*(oi-1)+2} = [idxs{oi} ' $\tilde{A}$\textbackslash$b$'];
        hold on;
    end
    xlabel('k',fss,fs); legend(mylegds, 'interpreter','latex');
	ylabel(['1-norm error on ' num2str(v.avm) ' BC'],fss,fs); set(gca,fss,fs);
    
    figure;
    for oi = 1:nrObsts
        loglog(v.ks, v.compresErr(length(v.ks)*(oi-1)+(1:length(v.ks)),2), ll{oi}, lws, lw);
        hold on;
    end
    xlabel('k',fss,fs); legend(idxs, 'interpreter','latex');
    ylabel('$||\tilde{A}c-b||_2/||b||_2$','interpreter','latex',fss,fs); set(gca,fss,fs);
    
    figure;
    for oi = 1:nrObsts
        loglog(v.ks, v.perc(length(v.ks)*(oi-1)+(1:length(v.ks)),2), ll{oi}, lws, lw);
        hold on;
    end
    xlabel('k',fss,fs); legend(idxs, 'interpreter','latex');
	ylabel('% nonzeros',fss,fs); set(gca,fss,fs);
    
    mylegds = cell(2*nrObsts,1);
    figure;
    for oi = 1:nrObsts
        loglog(v.ks, v.conds(length(v.ks)*(oi-1)+(1:length(v.ks)),:), ll{oi}, lws, lw);
        mylegds{2*(oi-1)+1} = [idxs{oi} ' $A$'];
        mylegds{2*(oi-1)+2} = [idxs{oi} ' $\tilde{A}$'];
        hold on;
    end
    xlabel('k',fss,fs); legend(mylegds, 'interpreter','latex');
	ylabel('Estimated condition number',fss,fs); set(gca,fss,fs);
    return
end

% Plot the results for the given indices and one scattering obstacle.
if isfield(v, 'conds')
    figure;
    loglog(v.ks, v.conds(idxs,1), 'b+', mss, ms); hold on;
    loglog(v.ks, v.conds(idxs,2), 'rx', mss, ms); 
    xlabel('k',fss,fs); legend({'$A$', '$\tilde{A}$'}, 'interpreter','latex');
	ylabel('Condition number',fss,fs); set(gca,fss,fs);
end
if isfield(v, 'errBCavm')
    figure;
    for qwer = 1:size(v.errBCavm,2)
        loglog(v.ks, v.errBCavm(idxs,qwer), [l{qwer} c(qwer)], mss,ms); hold on;
    end
    xlabel('k',fss,fs); 
	legend({'$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$','GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
        'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
		'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
	ylabel(['1-norm error on ' num2str(v.avm) ' BC'],fss,fs); set(gca,fss,fs);
end
if isfield(v, 'errTrueF')
	figure; loglog(v.ks, v.errTrueF(idxs,:), lws, lw); xlabel('k',fss,fs); 
	legend({'$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$','GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
        'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
		'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
	ylabel('Error wrt analytic solution',fss,fs); set(gca,fss,fs);
end

if isfield(v, 'nnz')
	figure;
    loglog(v.ks, v.nnz(idxs,1), 'b-.', lws,lw); hold on;
    loglog(v.ks, v.nnz(idxs,2), 'r*', mss,ms); hold on;
	loglog(v.ks, v.nnz(idxs(end),2)*(v.ks/v.ks(end)).^(3/2), 'g--', lws, lw); xlabel('k',fss,fs); 
	legend({'\# elements', '\# nonzeros','$\mathcal{O}(k^{3/2})$'}, 'interpreter', 'latex'); set(gca,fss,fs);
end
if isfield(v, 'perc')
	figure; loglog(v.ks, v.perc(idxs,:), lws, lw); hold on;
	loglog(v.ks, v.perc(idxs(end),2)*(v.ks/v.ks(end)).^(-1/2), 'r--', lws, lw); xlabel('k',fss,fs); 
	legend({'\% elements', '\% nonzeros','$\mathcal{O}(k^{-1/2})$'}, 'interpreter', 'latex'); set(gca,fss,fs);
end
if isfield(v,'errSol')
	figure; 
    for qwer = 2:size(v.errSol,2)
        loglog(v.ks, v.errSol(idxs,qwer), [l{qwer} c(qwer)], mss,ms); hold on;
    end
    xlabel('k','FontSize',fs);
	legend({'$\tilde{A}$\textbackslash$b$','GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
        'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
		'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
	ylabel('$|| \tilde{c}-$c$ ||/|| $c$ ||$','interpreter','latex',fss,fs); set(gca,fss,fs);
end
if isfield(v, 'errBCcol')
	figure; loglog(v.ks, v.errBCcol(idxs,2:end), lws, lw); xlabel('k',fss,fs); 
	legend({'$\tilde{A}$\textbackslash$b$','GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
        'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
		'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
	ylabel('Error on BC on collocation points',fss,fs); set(gca,fss,fs);
end
if isfield(v, 'compresErr')
	figure; loglog(v.ks, v.compresErr(idxs,2:end), lws, lw); xlabel('k',fss,fs); 
	ylabel('$||\tilde{A}c-b||_2/||b||_2$','interpreter','latex',fss,fs); set(gca,fss,fs);
end

if isfield(v, 'nbIter')
    figure;  
    for qwer = 1:size(v.nbIter,2)
        loglog(v.ks,v.nbIter(idxs,qwer), [l{qwer} c(qwer)], mss, ms); hold on;
    end
    xlabel('k',fss,fs); set(gca,fss,fs);
	legend({'$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$','$(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
        '$A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs); ylabel('# Gmres iterations',fss,fs);
end
if isfield(v,'timeA') && isfield(v,'timeSol')
    figure;
    if (size(v.timeA,2) >=4)
        lw = 2.3;
        loglog(v.ks, v.timeA(idxs,1:2), lws, lw); 
        xlabel('k',fss,fs); 
        ylabel('Time(s)');
        set(gca,fss,fs);
        nzTiCo = find(v.timeA(idxs,4));
        hold on; 
        loglog(v.ks, v.timeSol(idxs,1:3), '-.', lws, lw);
        loglog(v.ks, v.timeSol(idxs,4:end), '--', lws, lw);
        loglog(v.ks(nzTiCo), v.timeA(idxs(nzTiCo),4), 'k*', 'MarkerSize', 12);
        legend({'$A$','$\tilde{A}$','$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$',...
            'GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$', '$R_{m,n}$'}, 'interpreter','latex',fss,fs);
        
        % Above figure for article
        figure;
        loglog(v.ks, v.timeA(idxs,1:3), lws, lw); xlabel('k',fss,fs); set(gca,fss,fs);
        ylabel('Time(s)');
        nzTiCo = find(v.timeA(idxs,4));
        hold on; loglog(v.ks(nzTiCo), v.timeA(idxs(nzTiCo),4), 'k*', 'MarkerSize', 12);
        loglog(v.ks, v.timeSol(idxs,1:3), '-.', lws, lw);
        loglog(v.ks, v.timeSol(idxs,4:end), '--', lws, lw);
        legend({'$A$','$\tilde{A}$', 'Total','Correlations','$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$',...
            'GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
            'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
            'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
    else
        loglog(v.ks, v.timeA(idxs,:), lws, lw); set(gca,fss,fs);
        if size(v.timeSol,2) > 3
            loglog(v.ks, v.timeSol(idxs,1:3), '-.', lws, lw);
            loglog(v.ks, v.timeSol(idxs,4:6), '--', lws, lw);
        else
            loglog(v.ks, v.timeSol(idxs,:), lws, lw);
        end
        legend({'$A$','$\tilde{A}$', 'Total','$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$',...
            'GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
            'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
            'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs);
    end
    ylabel('Computation time (s)','FontSize',fs);
    xlabel('k',fss,fs); set(gca,fss,fs);
elseif isfield(v,'timeA')
    if (size(v.timeA,2) >=4)
        figure; loglog(v.ks, v.timeA(idxs,1:3), lws, lw); xlabel('k',fss,fs); set(gca,fss,fs);
        nzTiCo = find(v.timeA(idxs,4));
        hold on; loglog(v.ks(nzTiCo), v.timeA(idxs(nzTiCo),4), 'k*', 'MarkerSize', 12);
    else
        figure; loglog(v.ks, v.timeA(idxs,:), lws, lw); xlabel('k',fss,fs); set(gca,fss,fs);
    end
	legend({'$A$','$\tilde{A}$', 'Total','Correlations'}, 'interpreter','latex'); ylabel('Computation time (s)','FontSize',fs);
elseif isfield(v,'timeSol')
	figure; loglog(v.ks, v.timeSol(idxs,:), lws, lw); xlabel('k',fss,fs); set(gca,fss,fs);
	legend({'$A$\textbackslash$b$','$\tilde{A}$\textbackslash$b$','GM $A$\textbackslash$b$','GM $\tilde{A}$\textbackslash$b$',...
        'GM $(\tilde{A}^{-1}A)$\textbackslash$(\tilde{A}^{-1}b)$',...
		'GM $A$\textbackslash$b$ with $x_0=\tilde{x}$'}, 'interpreter','latex',fss,fs); ylabel('Time for c (s)',fss,fs);
end
if isfield(v,'field')
	figure;
    pcolor(v.xs, v.ys, abs(v.field))
    shading interp;
	h(3) = colorbar;
	ylabel(h(3),'Absolute value of the total solution')
	hold on;
	if iscell(v.parametr)
		for obst = 1:length(v.parametr)
			vp = v.parametr{obst};
			plot3(vp(1,:), vp(2,:), max(max(abs(v.field)))*ones(size(vp(1,:))), 'w', lws, lw);
		end
	else
		plot3(v.parametr(1,:), v.parametr(2,:), max(max(abs(v.field)))*ones(size(v.parametr(1,:))), 'w', lws, lw);
    end
    set(gca,fss,fs);
end
if exist('A','var') && ~isempty(A)
	figure; spy(A);
	ylabel('Collocation \tau', fss, fs);
	xlabel('Interaction \tau', fss, fs);
	set(gca,fss, fs)
	set(gca,'XTick',0:round(size(A,1)/5):size(A,1))
	set(gca,'XTickLabel',0:0.2:1)
	set(gca,'YTick',0:round(size(A,2)/10):size(A,2))
	set(gca,'YTickLabel',0:0.1:1)
    set(gca,'tickdir','out')
    set(gca, 'Layer','top')
end
