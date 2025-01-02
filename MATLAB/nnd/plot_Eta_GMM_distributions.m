function plot_Eta_GMM_distributions(vEta,MM_models,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 2 June 2022
%   ...
%   version 1.0.1, 29 August 2024
%
    sFitModel = 'GMM'; % 'GMM', 'Weibul' 
    nBins     = 30;
    sBinType  = 'loglog';
    bSaveFig  = true;
    for k = 1:length(varargin)
        if strcmp('FitModel',varargin{k})
            sFitModel = varargin{k+1};
        end
        if strcmp('Binning',varargin{k})
            sBinType = varargin{k+1};
        end
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end
    sPlateNum = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)'};

    nEtabins = 40;
    yl2 = 0.7;
    sLsty = {'--','-.',':'};
    % compute the histogram
    %vEta_pdf = histcounts(log10(eta),vEedges,'Normalization','pdf');        % histogram of the rescaled proximity
    [vEta_pdf, vEedges] = histcounts(log10(vEta),nEtabins,'Normalization','pdf'); % histogram of the rescaled proximity

    fwx = 900;
    fwy = 800;
    if Model.nGMMmax == 2
        ny = 1;
        nx = 2;
        fwy = 450;
    elseif Model.nGMMmax == 3 || Model.nGMMmax == 4
        ny = 2;
        nx = 2;
    elseif Model.nGMMmax == 5 || nPar == 6
        ny = 2;
        nx = 3;
    elseif Model.nGMMmax > 6
        nx = 3;
        ny = 3;
    end
    
    figure('Name','Distribution of the rescaled distances','Position',[550 300 fwx fwy]);
    set(gcf,'Color','w'); % background color for the figure
    for nc = 1:Model.nGMMmax
        subplot(ny,nx,nc)
        % plot the \eta distribution
        %histogram(log10(eta),nEtabins,'Normalization','pdf');
        histogram('BinEdges',vEedges,'BinCounts',vEta_pdf);
        ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
        ax.YLim = [0, yl2];
        grid on
        hold on
        %axis tight
        axis manual
        vXmodel = linspace(vEedges(1),vEedges(end),200)';
        cLegendEta{1} = 'Density of log_{10}(\eta)';
        if strcmp(sFitModel,'GMM')
            cLegend = plot_gmm_fit(MM_models{nc},vXmodel,'Color','#A2142F','Components');
        elseif strcmp(sFitModel,'Weibul')
        end
        cLegendEta(2:length(cLegend)+1) = cLegend;
        legend(cLegendEta,'Location','northwest','FontSize',7);
        legend('boxoff');
        xlabel('Nearest-neighbour proximity, $\log_{10}(\eta)$','Interpreter','latex');
        ylabel('Density','Interpreter','latex');
        text(0.9,0.95,sPlateNum{nc},'FontSize',12,'Units','normalized');
        hold off
    end
    sgtitle(Model.sTitle,'FontSize',12);

    if bSaveFig
        sName = [Model.sFileName,'_Eta_GMMs'];
        save_cf(gcf,sName,'fig','png','pdf');
    end
end
