function ax = plot_TR_declustered(T,R,BkgrdInds,AshkInds,crossBkgrdInds,crossAshkInds,NNDpar,options)
%
%   Plots the distribution of R vs T for clustered and background events
%   T, R - rescaled time and distance
%   BkgrdInds, AshkInds - indeces for background and aftershock events 
%   Model - structure that defines the seismicity model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 14 October 2024
%   ...
%   version 1.2.0, 27 October 2024
%
    arguments
        T double
        R double
        BkgrdInds logical
        AshkInds logical
        crossBkgrdInds logical
        crossAshkInds logical
        NNDpar struct
        options.Synthetic logical = false   % wether synthetic or real dataset 
        options.Colors double = []
        options.MarkerSize double = 10
        options.Title char = []          % to plot the title or not 
        options.Legend char = 'on'       % 'on'/'off' to plot the legend or not 
        options.YAxisLabel char = 'on'   % 'on'/'off' to plot the Y axis label or not 
        options.ThresholdValue char = 'on' % 'on'/'off' to plot the \eta value or not 
    end

    % Set plotting visuals
    if isempty(options.Colors)
        clrs = lines(5);
        colors = [clrs(1,:);clrs(2,:);clrs(1,:);clrs(2,:);[0 0 0]];
    else
        colors = options.Colors;
    end
    MkrSize = options.MarkerSize;
    if options.Synthetic
        crossMkrSize = MkrSize + 10;
        MkrAlpha = 0.2;
        crossMkrAlpha = 0.9;
        MrkEdgeCol = 'none';
        MrkEdgeAlpha = 0.0;
    else
        crossMkrSize = MkrSize + 20;
        MkrAlpha = 0.6;
        crossMkrAlpha = 0.8;
        MrkEdgeCol = 'k';
        MrkEdgeAlpha = 0.3;
    end
    txtSize = 11;
    LabelFontSize = 12;
    LgdFontSize = 12;
    TitleFontSize = 14;
    TRhistXLim = [log10(NNDpar.TLim), 2];
    TRhistYLim = [log10(NNDpar.RLim), 2];
    
    % Plot NND distribution
    p(1) = scatter(T(BkgrdInds),R(BkgrdInds),MkrSize,colors(1,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Background events');
    hold on
    p(2) = scatter(T(AshkInds),R(AshkInds),MkrSize,colors(2,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Aftershocks');
    if ~isempty(crossBkgrdInds)
        scatter(T(crossBkgrdInds),R(crossBkgrdInds),crossMkrSize,colors(3,:),'p','filled',...
            'MarkerEdgeColor',MrkEdgeCol,'MarkerEdgeAlpha',MrkEdgeAlpha,'MarkerFaceAlpha',crossMkrAlpha,'HandleVisibility','off');
    end
    if ~isempty(crossAshkInds)
        scatter(T(crossAshkInds),R(crossAshkInds),crossMkrSize,colors(4,:),'p','filled',...
            'MarkerEdgeColor',MrkEdgeCol,'MarkerEdgeAlpha',MrkEdgeAlpha,'MarkerFaceAlpha',crossMkrAlpha,'HandleVisibility','off');
    end
    p(3) = plot(TRhistXLim(1:2),NNDpar.vLog10EtaThresh-TRhistXLim(1:2),'--','Color',colors(5,:),'LineWidth',1.5,'DisplayName','NND threshold');
    hold off
    box on
    grid on
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.XLim = TRhistXLim(1:2);
    ax.YLim = TRhistYLim(1:2);
    %ax.XTick = round(ax.XLim(1),0):TRhistXLim(3):round(ax.XLim(2),0);
    %ax.YTick = round(ax.YLim(1),0):TRhistYLim(3):round(ax.YLim(2),0);
    if strcmp(options.ThresholdValue,'on')
        text(0.05,0.05,['NND threshold: $\log_{10}(\eta) = ',num2str(NNDpar.vLog10EtaThresh),'$'],'Units','normalized','Interpreter','latex','FontSize',LgdFontSize);
    end
    if strcmp(options.Legend,'on')
        legend(p,'Location','NW','FontSize',LgdFontSize,'Interpreter','latex');
    end
    xlabel('Rescaled time, $\log_{10}(T)$','FontSize',LabelFontSize,'Interpreter','latex');
    if strcmp(options.YAxisLabel,'on')
        ylabel('Rescaled distance, $\log_{10}(R)$','FontSize',LabelFontSize,'Interpreter','latex');
    end
    title(options.Title,'FontSize',TitleFontSize,'Interpreter','latex');
end

