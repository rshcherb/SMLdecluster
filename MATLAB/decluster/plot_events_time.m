function ax = plot_events_time(vCat,BkgrdInds,AshkInds,Model,OptArgs)
%
%   Plots the evolution in time of the sequence for clustered and background events
%   vCat - catalogue of events
%   BkgrdInds, AshkInds - indeces for background and aftershock events 
%   Model - structure that defines the seismicity model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 16 October 2024
%   ...
%   version 1.1.0, 3 December 2024
%
    arguments
        vCat double
        BkgrdInds logical
        AshkInds logical
        Model struct
        OptArgs.Title char = []   % to plot the title or not 
        OptArgs.Legend char = 'on'   % 'on'/'off' to plot the legend or not 
        OptArgs.YAxisLabel char = 'on'   % 'on'/'off' to plot the Y axis label or not 
        OptArgs.DateFormat char = 'date'   % 'date', 'day' format to use for x-axis 
        OptArgs.FigLabel cell = {'(a)', '(b)'} % letter lables for each figure
    end

    BkgrdLat = vCat(BkgrdInds,2);
    AshkLat = vCat(AshkInds,2);
    vTbgrd = vCat(BkgrdInds,1);
    vTashk = vCat(AshkInds,1);
    if strcmp(OptArgs.DateFormat,'date')
        vTbgrd = datenum(Model.vDateStart) + vTbgrd;
        vTashk = datenum(Model.vDateStart) + vTashk;
    end
    BkgrdSize = (Model.fMc - 1.0 + vCat(BkgrdInds,4) - min(vCat(BkgrdInds,4))).^3;
    AshkSize = (Model.fMc - 1.0 + vCat(AshkInds,4) - min(vCat(AshkInds,4))).^3;

    % Set plotting visuals
    clrs = lines(5);
    colors = [clrs(1,:);clrs(2,:);[0 0 0]];
    MkrAlpha = 0.6;
    LabelFontSize = 12;
    LgdFontSize = 10;
    TitleFontSize = 14;
    BkgrdAlpha = 1.0;
    AshkAlpha = 1.0;
    xlbl = -0.06; ylbl = 1.05;

    tiledlayout(2,1,TileSpacing="tight");

    ax(1) = nexttile;

    p1(1) = scatter(vTbgrd,BkgrdLat,BkgrdSize,colors(1,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Background events');
    hold on
    p1(2) = scatter(vTashk,AshkLat,AshkSize,colors(2,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Aftershocks');
    ax(1) = gca; ax(1).XMinorTick = 'on'; ax(1).YMinorTick = 'on';
    hold off
    grid on
    box on
    axis tight
    text(xlbl,ylbl,OptArgs.FigLabel{1},'FontSize',15,'Units','normalized','Interpreter','latex');
    
    if strcmp(OptArgs.Legend,'on')
        legend(p1,'Location','NW','FontSize',LgdFontSize,'Interpreter','latex');
    end
    if strcmp(OptArgs.DateFormat,'date')
        xlabel('Date (year)','FontSize',LabelFontSize,'Interpreter','latex');
        datetick('x','yyyy','keeplimits');
    elseif strcmp(OptArgs.DateFormat,'day')
        xlabel('Time (days)','FontSize',LabelFontSize,'Interpreter','latex');
    end
    if strcmp(OptArgs.YAxisLabel,'on')
        ylabel('Latitude','FontSize',LabelFontSize,'Interpreter','latex');
    end
    title(OptArgs.Title,'FontSize',TitleFontSize,'Interpreter','latex');

    ax(2) = nexttile;

    p2(1) = scatter(vTbgrd,BkgrdLat,BkgrdSize,colors(1,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Background events');
    hold on
    %p2(2) = scatter(vTashk,AshkLat,AshkSize,colors(2,:),'filled','MarkerFaceAlpha',MkrAlpha,'DisplayName','Aftershocks');
    ax(2) = gca; ax(2).XMinorTick = 'on'; ax(2).YMinorTick = 'on';
    hold off
    grid on
    box on
    axis tight
    text(xlbl,ylbl,OptArgs.FigLabel{2},'FontSize',15,'Units','normalized','Interpreter','latex');
    
    if strcmp(OptArgs.Legend,'on')
        legend(p2,'Location','NW','FontSize',LgdFontSize,'Interpreter','latex');
    end
    if strcmp(OptArgs.DateFormat,'date')
        xlabel('Date (year)','FontSize',LabelFontSize,'Interpreter','latex');
        datetick('x','yyyy','keeplimits');
    elseif strcmp(OptArgs.DateFormat,'day')
        xlabel('Time (days)','FontSize',LabelFontSize,'Interpreter','latex');
    end
    if strcmp(OptArgs.YAxisLabel,'on')
        ylabel('Latitude','FontSize',LabelFontSize,'Interpreter','latex');
    end
end

