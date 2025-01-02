function plot_cluster_stats(NNDStats,vCat,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 3 June 2022
%   ...
%   version 1.2.1, 29 August 2024
%
    nBins       = [40, 40];
    bSaveFig    = true;
    for k = 1:length(varargin)
        if strcmp('NumBins',varargin{k})
            nBins = varargin{k+1};
        end
    end
    sPlateNum = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)'};

    nComp_indx = NNDStats.nComp_indx;
    nComp_size = NNDStats.nComp_size;
    ClstProp   = NNDStats.ClstProp;
    nClst      = NNDStats.nClst;      % number of clusters
    nSingle    = NNDStats.nSingle;    % number of singleton events

    % 
    figure('Name','Cluster component numbers','Position',[550 300 900 400]);
    set(gcf,'Color','w'); % background color for the figure

    subplot(1,2,1);
    bar(nComp_size,'EdgeColor','#0072BD','LineWidth',1.2);
    ax = gca; ax.YMinorTick = 'on';
    hold on
    xlabel('Cluster component #');
    ylabel('Number of events');
    text(0.9,0.95,'a)','FontSize',12,'Units','normalized');
    hold off

    subplot(1,2,2);
    %histogram(nComp_size,'BinMethod','integers');
    histogram(NNDStats.cClstComps);
    ax = gca; ax.YMinorTick = 'on';
    ax.YScale = 'log';
    ax.YLim = [0.7 Inf];
    axis padded
    xtickangle(ax,0);
    hold on
    xlabel('Number of events');
    ylabel('Number of family trees');
    text(0.9,0.95,'b)','FontSize',12,'Units','normalized');
    hold off
    sgtitle(Model.sTitle,'FontSize',12);

    if bSaveFig
        sName = [Model.sFileName,'_family_tree_stats'];
        save_cf(gcf,sName,'fig','png','pdf');
    end
end

