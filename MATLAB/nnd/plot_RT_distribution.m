function plot_RT_distribution(vT,vR,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.1.0, 2 December 2022
%
    sFitModel = 'GMM'; % 'GMM', 'Weibul' 
    nBins     = [40, 40];
    bSaveFig  = true;
    for k = 1:length(varargin)
        if strcmp('FitModel',varargin{k})
            sFitModel = varargin{k+1};
        end
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end

    sLsty = {'--','-.',':'};
    % compute the histogram
    %[mRT, Tedges, Redges] = histcounts2(log10(vT),log10(vR),nBins,'Normalization','pdf');
    [mRT, Tedges, Redges] = histcounts2(log10(vT),log10(vR),nBins);
    vXc = Tedges(1:end-1);
    vYc = Redges(1:end-1);
    %disp(mRT)
    mRT = mRT';
%
    figure('Name','2D distributions of the rescaled distances','Position',[550 300 700 600]);
    set(gcf,'Color','w'); % background color for the figure

    pcolor(vXc,vYc,log10(mRT));
    %pcolor(vXc,vYc,mRT);
    shading flat;
    %shading interp;
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    hold on
    axis manual
    colorbar;
    vXthr = vXc;
    cLegend{1} = '';
    for k = 1:length(Model.EtaThresh)
        vYthr = log10(Model.EtaThresh(k)) - vXthr;
        plot(vXthr,vYthr,'w','LineWidth',1.5,'LineStyle',sLsty{k});
        cLegend{k+1} = ['log_{10}(\eta_{thresh}) = ',num2str(log10(Model.EtaThresh(k)))];
    end
    legend(cLegend,'Location','southwest','FontSize',10,'TextColor','w');
    legend('boxoff');
    xlabel('Rescaled time, $\log_{10}(T)$','Interpreter','latex');
    ylabel('Rescaled distance, $\log_{10}(R)$','Interpreter','latex');
    title(Model.sTitle);
    hold off

    if bSaveFig
        save_cf(gcf,[Model.sFileName,'_TRhist_pcolor'],'fig','png','pdf');
    end
%
    figure('Name','2D distributions of the rescaled distances histogram','Position',[550 300 700 600]);
    set(gcf,'Color','w'); % background color for the figure

    histogram2(log10(vT),log10(vR),nBins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none','Normalization','pdf');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    xlim = ax.XLim;
    ylim = ax.YLim;
    hold on
    axis manual
    colorbar;
    vXthr = vXc;
    cLegend{1} = '';
    for k = 1:length(Model.EtaThresh)
        vYthr = log10(Model.EtaThresh(k)) - vXthr;
        plot(vXthr,vYthr,'w','LineWidth',1.5,'LineStyle',sLsty{k});
        cLegend{k+1} = ['log_{10}(\eta_{thresh}) = ',num2str(log10(Model.EtaThresh(k)))];
    end
    legend(cLegend,'Location','southwest','FontSize',10,'TextColor','w');
    legend('boxoff');
    xlabel('Rescaled time, $\log_{10}(T)$','Interpreter','latex');
    ylabel('Rescaled distance, $\log_{10}(R)$','Interpreter','latex');
    title(Model.sTitle);
    hold off

    if bSaveFig
        save_cf(gcf,[Model.sFileName,'_TRhist_histogram2'],'fig','png','pdf');
    end
%
    xgrid = linspace(xlim(1),xlim(2),nBins(1));
    ygrid = linspace(ylim(1),ylim(2),nBins(2));
    [x1,y1] = meshgrid(xgrid, ygrid);
    xi = [x1(:) y1(:)];
    [df, xep] = ksdensity([log10(vT), log10(vR)],xi);
    X = reshape(xep(:,1),length(xgrid),length(ygrid));
    Y = reshape(xep(:,2),length(xgrid),length(ygrid));
    Z = reshape(df,length(xgrid),length(ygrid));

    figure('Name','2D distributions of the rescaled distances: ksdensity','Position',[550 300 700 600]);
    set(gcf,'Color','w'); % background color for the figure

    contourf(X,Y,Z,15,'LineStyle','none');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    hold on
    axis manual
    plot(log10(vT),log10(vR),'.','Color','black');
    colorbar;
    vXthr = vXc;
    cLegend = {};
    for k = 1:length(Model.EtaThresh)
        vYthr = log10(Model.EtaThresh(k)) - vXthr;
        vgpl(k) = plot(vXthr,vYthr,'w','LineWidth',1.5,'LineStyle',sLsty{k});
        cLegend{k} = ['log_{10}(\eta_{thresh}) = ',num2str(log10(Model.EtaThresh(k)))];
    end
    legend(vgpl,cLegend,'Location','southwest','FontSize',10,'TextColor','w');
    legend('boxoff');
    xlabel('Rescaled time, $\log_{10}(T)$','Interpreter','latex');
    ylabel('Rescaled distance, $\log_{10}(R)$','Interpreter','latex');
    title(Model.sTitle);
    hold off

    if bSaveFig
        save_cf(gcf,[Model.sFileName,'_TRdensity_ksdensity'],'fig','png','pdf');
    end
end

