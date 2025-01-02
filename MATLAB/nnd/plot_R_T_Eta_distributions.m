function plot_R_T_Eta_distributions(vEta,vT,vR,MM_Eta,MM_T,MM_R,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.1.0, 17 October 2024
%
    nBins     = 30;
    sBinType  = 'loglog';
    bSaveFig  = true;
    for k = 1:length(varargin)
        if strcmp('Binning',varargin{k})
            sBinType = varargin{k+1};
        end
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end

    nRbins = 30;
    nTbins = 30;
    nEtabins = 70;
    yl2 = 0.65;
    sLsty = {'--','-.',':'};
    % compute the histogram
    %vR_pdf = histcounts(log10(vR),vRedges,'Normalization','pdf');           % histogram of the rescled distances
    [vR_pdf, vRedges] = histcounts(log10(vR),nRbins,'Normalization','pdf');      % histogram of the rescaled distances
    %vT_pdf = histcounts(log10(vT),vTedges,'Normalization','pdf');           % histogram of the rescaled times
    [vT_pdf, vTedges] = histcounts(log10(vT),nTbins,'Normalization','pdf');      % histogram of the rescaled times
    %vEta_pdf = histcounts(log10(eta),vEedges,'Normalization','pdf');        % histogram of the rescaled proximity
    [vEta_pdf, vEedges] = histcounts(log10(vEta),nEtabins,'Normalization','pdf'); % histogram of the rescaled proximity

    figure('Name','Distributions of the rescaled distances','Position',[550 300 900 800]);
    set(gcf,'Color','w'); % background color for the figure

    % plot the R distribution
    subplot(2,2,1)
%         if strcmp(sBinType,'linear')
%             semilogy(vX,vDt_pdf,'s','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',8);
%         elseif strcmp(sBinType,'loglog')
%             loglog(vX,vDt_pdf,'s','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120','MarkerSize',8);
%         end
    %histogram(log10(vR),nRbins,'Normalization','pdf');
    histogram('BinEdges',vRedges,'BinCounts',vR_pdf);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.YLim = [0,yl2];
    grid on
    hold on
    %axis manual
    %axis tight
    cLegendR{1} = 'Density of log_{10}(R)';
    vXmodel = linspace(vRedges(1),vRedges(end),200)';
    if strcmp(Model.NNDpar.MixtModel,'GMM')
        cLegend = plot_gmm_fit(MM_R,vXmodel,'Color','#A2142F','NumCompMax',Model.nGMMmax);
    elseif strcmp(Model.NNDpar.MixtModel,'WMM') % Weibull mixture model
    end
    cLegendR(2:length(cLegend)+1) = cLegend;
    legend(cLegendR,'Location','northwest','FontSize',6);
    legend('boxoff');
    xlabel('Rescaled distance, $\log_{10}(R)$','Interpreter','latex');
    ylabel('Density','Interpreter','latex');
    text(0.91,0.95,'a)','FontSize',12,'Units','normalized');
    hold off

    % plot the T distribution
    subplot(2,2,2)

    %histogram(log10(vT),nTbins,'Normalization','pdf');
    histogram('BinEdges',vTedges,'BinCounts',vT_pdf);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    grid on
    hold on
    %axis manual
    %axis tight
    vXmodel = linspace(vTedges(1),vTedges(end),200)';
    if strcmp(Model.NNDpar.MixtModel,'GMM')
        cLegend = plot_gmm_fit(MM_T,vXmodel,'Color','#A2142F');
    elseif strcmp(Model.NNDpar.MixtModel,'Weibul')
    end
    cLegendT{1} = 'Density of log_{10}(T)';
    cLegendT(2:length(cLegend)+1) = cLegend;
    legend(cLegendT,'Location','northwest','FontSize',6);
    legend('boxoff');
    xlabel('Rescaled time, $\log_{10}(T)$','Interpreter','latex');
    ylabel('Density','Interpreter','latex');
    text(0.91,0.95,'b)','FontSize',12,'Units','normalized');
    hold off
    
    % plot the \eta distribution
    subplot(2,2,[3,4])

    %histogram(log10(eta),nEtabins,'Normalization','pdf');
    histogram('BinEdges',vEedges,'BinCounts',vEta_pdf);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    grid on
    hold on
    %axis tight
    axis manual
    vXmodel = linspace(vEedges(1),vEedges(end),200)';
    if strcmp(Model.NNDpar.MixtModel,'GMM')
        cLegend = plot_gmm_fit(MM_Eta,vXmodel,'Color','#A2142F','Components');
    elseif strcmp(Model.NNDpar.MixtModel,'Weibul')
    end
    cLegendEta{1} = 'Density of log_{10}(\eta)';
    cLegendEta(2:length(cLegend)+1) = cLegend;
    nl = length(cLegendEta);
    for k = 1:length(Model.EtaThresh)
        plot([log10(Model.EtaThresh(k)), log10(Model.EtaThresh(k))],[ax.YLim(1), ax.YLim(2)],'k','LineWidth',1.5,'LineStyle',sLsty{k});
        nl = nl + 1;
        cLegendEta{nl} = ['log_{10}(\eta_{thresh}) = ',num2str(log10(Model.EtaThresh(k)))];
    end
    legend(cLegendEta,'Location','northwest','FontSize',7);
    legend('boxoff');
    xlabel('Nearest-neighbour proximity, $\log_{10}(\eta)$','Interpreter','latex');
    ylabel('Density','Interpreter','latex');
    hold off
    text(0.95,0.95,'c)','FontSize',12,'Units','normalized');
    sgtitle(Model.sTitle,'FontSize',12);

    if bSaveFig
        sName = [Model.sFileName,'_TReta'];
        save_cf(gcf,sName,'fig','png','pdf');
    end
end
