function plot_NND_T_R(testCats,realCat,NNDpar)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        testCats table
        realCat table
        NNDpar struct
    end

    % Set plotting visuals
    colors = [rgb('black');rgb('gray');[0.3010 0.7450 0.9330]];
    realCatMkrAlpha = 1;
    testCatsMkrAlpha = 0.2;
    MkrSize = 5;
    TRhistXLim = [log10(NNDpar.TLim), 2];
    TRhistYLim = [log10(NNDpar.RLim), 2];
    LabelFontSize = 14;
    LgdFontSize = 12;

    scatter(testCats.T1,testCats.R1,MkrSize,colors(2,:),'filled','MarkerFaceAlpha',testCatsMkrAlpha,'DisplayName','Synthetic Events');
    hold on
    scatter(realCat.T1,realCat.R1,MkrSize,colors(1,:),'filled','MarkerFaceAlpha',realCatMkrAlpha,'DisplayName','Real Events');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';

    nBins = [80, 80];
    xlim = ax.XLim;
    ylim = ax.YLim;
    xgrid = linspace(xlim(1),xlim(2),nBins(1));
    ygrid = linspace(ylim(1),ylim(2),nBins(2));
    [x1,y1] = meshgrid(xgrid, ygrid);
    xi = [x1(:) y1(:)];
    [df, xep] = ksdensity([testCats.T1, testCats.R1],xi);
    X = reshape(xep(:,1),length(xgrid),length(ygrid));
    Y = reshape(xep(:,2),length(xgrid),length(ygrid));
    Z = reshape(df,length(xgrid),length(ygrid));
    contour(X,Y,Z,7,'LineColor','w','LineWidth',0.5,'LineStyle','-','HandleVisibility','Off');
%     [h,XEdges,YEdges] = histcounts2(testCats.T1,testCats.R1,40);
%     xvals = XEdges(1:end-1) + ((XEdges(2)-XEdges(1))/2);
%     yvals = YEdges(1:end-1) + ((YEdges(2)-YEdges(1))/2);
%     contour(xvals,yvals,flipud(rot90(h)),5,'LineColor','w','LineWidth',0.1,'LineStyle','-','HandleVisibility','Off');
    plot(TRhistXLim(1:2),NNDpar.vLog10EtaThresh-TRhistXLim(1:2),'--','Color',colors(3,:),'LineWidth',1.5,'DisplayName','NND Threshold');
    hold off
    box on
    grid on
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.XLim = TRhistXLim(1:2);
    ax.YLim = TRhistYLim(1:2);
    ax.XTick = round(ax.XLim(1),0):TRhistXLim(3):round(ax.XLim(2),0);
    ax.YTick = round(ax.YLim(1),0):TRhistYLim(3):round(ax.YLim(2),0);
    xlabel('Rescaled time, $\log_{10}(T)$','FontSize',LabelFontSize,'Interpreter','latex')
    ylabel('Rescaled distance, $\log_{10}(R)$','FontSize',LabelFontSize,'Interpreter','latex')
    legend('Location','northeast','box','off','FontSize',LgdFontSize,'Interpreter','latex')
end