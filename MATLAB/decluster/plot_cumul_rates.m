function ax = plot_cumul_rates(vCat,BkgrdInds,AshkInds,vDateStart,options)
%
%   Plots the cumulative event rates for the clustered and background events
%   vCat - catalogue of events
%   BkgrdInds, AshkInds - indeces for background and aftershock events 
%   Model - structure that defines the seismicity model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 14 October 2024
%   ...
%   version 1.1.0, 27 October 2024
%
    arguments
        vCat double
        BkgrdInds logical
        AshkInds logical
        vDateStart double
        options.Title char = []   % to plot the title or not 
        options.Legend char = 'on'   % 'on'/'off' to plot the legend or not 
        options.YAxisLabel char = 'on'   % 'on'/'off' to plot the Y axis label or not 
        options.DateFormat char = 'date'   % 'date', 'day' format to use for x-axis 
        options.PieChart char = 'on'   % 'on'/'off' to plot the piechart or not 
        options.PieSize double = 0.1;
        options.TestResults struct = []   % 'on'/'off' to plot the test results or not 
    end

    BkgrdRate = cumsum(ones(sum(BkgrdInds),1))/sum(BkgrdInds);
    AshkRate = cumsum(ones(sum(AshkInds),1))/sum(AshkInds);
    vTbgrd = vCat(BkgrdInds,1);
    vTashk = vCat(AshkInds,1);
    t0 = vCat(find(BkgrdInds,1),1);
    t1 = vCat(find(BkgrdInds,1,'last'),1);
    if strcmp(options.DateFormat,'date')
        vTbgrd = datenum(vDateStart) + vTbgrd;
        vTashk = datenum(vDateStart) + vTashk;
        t0 = datenum(vDateStart) + t0;
        t1 = datenum(vDateStart) + t1;
    end

    % Set plotting visuals
    clrs = lines(5);
    colors = [clrs(1,:);clrs(2,:);[0 0 0]];
    txtSize = 11;
    LabelFontSize = 12;
    LgdFontSize = 12;
    TitleFontSize = 14;
    BkgrdAlpha = 1.0;
    AshkAlpha = 1.0;
    BkgrdLineWidth = 1.6;
    AshkLineWidth = 1.5;
    
    p(1) = plot(vTbgrd,BkgrdRate,'-','LineWidth',BkgrdLineWidth,'Color',[colors(1,:) BkgrdAlpha],'DisplayName','Background rate');
    hold on
    p(2) = plot(vTashk,AshkRate,'-','LineWidth',AshkLineWidth,'Color',[colors(2,:) AshkAlpha],'DisplayName','Aftershock rate');
    plot([t0,t1],[0,1],'k--','LineWidth',1.);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    if ~isempty(options.TestResults)
        % text(0.6,0.30,['Bridge test $p$-value: ',num2str(options.TestResults.Bridge.p)],'Units','normalized','Interpreter','latex','FontSize',LgdFontSize);
        text(0.6,0.25,['KS test $p$-value:     ',num2str(options.TestResults.KS.p)],'Units','normalized','Interpreter','latex','FontSize',LgdFontSize);
        text(0.6,0.20,['BZ test $p$-value:     ',num2str(options.TestResults.BZ.p)],'Units','normalized','Interpreter','latex','FontSize',LgdFontSize);
    end
    hold off
    grid on
    box on
    axis tight
    
    if strcmp(options.Legend,'on')
        legend(p,'Location','NW','FontSize',LgdFontSize,'Interpreter','latex');
    end
    if strcmp(options.DateFormat,'date')
        xlabel('Date (year)','FontSize',LabelFontSize,'Interpreter','latex');
        datetick('x','yyyy','keeplimits');
    elseif strcmp(options.DateFormat,'day')
        xlabel('Time (days)','FontSize',LabelFontSize,'Interpreter','latex');
        %xlabel('Time (days from $T_s$)','FontSize',LabelFontSize,'Interpreter','latex');
    end
    if strcmp(options.YAxisLabel,'on')
        ylabel('Cumulative rate','FontSize',LabelFontSize,'Interpreter','latex');
    end
    if strcmp(options.PieChart,'on')
        plot_bkgrd_ashk_piechart(sum(AshkInds),sum(BkgrdInds),'ax',ax,'PieSize',options.PieSize);
    end
    title(options.Title,'FontSize',TitleFontSize,'Interpreter','latex');
end

