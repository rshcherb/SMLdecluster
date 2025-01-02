function ax = plot_eta_declustered(eta,BkgrdInds,AshkInds,NNDpar,OptArgs)
%
%   Plots the distribution of the proximity \eta for clustered and background events
%   eta - proximity distance \eta
%   BkgrdInds, AshkInds - indeces for background and aftershock events 
%   Model - structure that defines the seismicity model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 15 October 2024
%   ...
%   version 1.1.0, 27 October 2024
%
    arguments
        eta double
        BkgrdInds logical
        AshkInds logical
        NNDpar struct
        OptArgs.PDF char = 'Normal'   % 'Normal', 'Weibull' pdf to fit to the distributions 
        OptArgs.Legend char = 'on'   % 'on'/'off' to plot the legend or not 
        OptArgs.YAxisLabel char = 'on'   % 'on'/'off' to plot the Y axis label or not 
        OptArgs.Title char = []   % to plot the title or not 
        OptArgs.PieChart char = 'on'   % 'on'/'off' to plot the piechart or not 
        OptArgs.PieSize double = 0.09;
    end

    % Set plotting visuals
    clrs = lines(5);
    colors = [clrs(1,:);clrs(2,:);[0 0 0]];
    LabelFontSize = 12;
    LgdFontSize = 12;
    TitleFontSize = 14;
    EtahistXLim = [log10(NNDpar.EtaLim) 2];
    EtahistYLim = [NNDpar.EtaYLim 0.1];
    
    eta_Bkgrd = eta(BkgrdInds);
    eta_Bkgrd = eta_Bkgrd(eta_Bkgrd > EtahistXLim(1) & eta_Bkgrd < EtahistXLim(2));
    eta_Ashk = eta(AshkInds);
    eta_Ashk = eta_Ashk(eta_Ashk > EtahistXLim(1) & eta_Ashk < EtahistXLim(2));
    % histograms for \eta
    % set histogram properties
    nbins = 35;
    binWidth = 0.22;

    hBkgrd = histogram(eta_Bkgrd,nbins,'BinWidth',binWidth,'Normalization','pdf',...
        'FaceAlpha',0.5,'DisplayName','Background events','EdgeColor',colors(1,:));
    hold on
    hAshk = histogram(eta_Ashk,nbins,'BinWidth',binWidth,'Normalization','pdf',...
        'FaceAlpha',0.5,'DisplayName','Aftershocks','EdgeColor',colors(2,:));

    pdBkgrd = fitdist(eta_Bkgrd - min(hBkgrd.BinEdges),OptArgs.PDF);
    pdAshk = fitdist(eta_Ashk - min(hAshk.BinEdges),OptArgs.PDF);
    
    xbkgr = linspace(min(eta_Bkgrd),max(eta_Bkgrd),200);
    xashk = linspace(min(eta_Ashk),max(eta_Ashk),200);

    p(1) = plot(xbkgr,pdf(pdBkgrd,xbkgr - min(hBkgrd.BinEdges)),'k-','LineWidth',1.5,'DisplayName',[OptArgs.PDF,' PDF']);
    plot(xashk,pdf(pdAshk,xashk - min(hAshk.BinEdges)),'k-','LineWidth',1.5,'HandleVisibility','off');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.XLim = EtahistXLim(1:2);
    ax.YLim = EtahistYLim(1:2);
    p(2) = plot([NNDpar.vLog10EtaThresh,NNDpar.vLog10EtaThresh],[0,ax.YLim(2)],'--','Color',colors(3,:),...
        'LineWidth',1.5,'DisplayName',['NND threshold: $\log_{10}(\eta) = ',num2str(NNDpar.vLog10EtaThresh),'$']);
    if strcmp(OptArgs.Legend,'on')
        legend([hBkgrd, hAshk, p],'Location','NW','FontSize',LgdFontSize,'Interpreter','latex');
    end
    hold off
    grid on
    xlabel('Nearest-neighbor proximity, $\log_{10}(\eta)$','FontSize',LabelFontSize,'Interpreter','latex');
    if strcmp(OptArgs.YAxisLabel,'on')
        ylabel('Density','FontSize',LabelFontSize,'Interpreter','latex');
    end
    if strcmp(OptArgs.PieChart,'on')
        pax = plot_bkgrd_ashk_piechart(sum(AshkInds),sum(BkgrdInds),'ax',ax,'PieSize',OptArgs.PieSize);
        pax.Position = [ax.Position(1)+0.05*ax.Position(3),ax.Position(2)+0.6*ax.Position(4),OptArgs.PieSize,OptArgs.PieSize];
    end
    title(OptArgs.Title,'FontSize',TitleFontSize,'Interpreter','latex');
end

