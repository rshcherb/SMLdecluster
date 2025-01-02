function plot_NND_distribution(vEta,vT,vR)
% Plots the 1D and 2D distributions of the NND
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    subplot(2,1,1)
    if numel(vEta)<5000
        nbins = 30;
    else
        nbins = 40;
    end
    log10vEta = log10(vEta);
    log10vT = log10(vT);
    log10vR = log10(vR);
    real_idx = ~isinf(log10vEta) & ~isnan(log10vEta) & ~isinf(log10vT) & ~isnan(log10vT)...
        & ~isinf(log10vR) & ~isnan(log10vR);
    log10vEta = log10vEta(real_idx);
    
    log10vT = log10vT(real_idx);
    
    log10vR = log10vR(real_idx);
    plot_NND2Dhist(log10vT,log10vR,nbins)
    ax1 = gca;
    axlim = axis(ax1);
    log10T = linspace(axlim(1),axlim(2),100);
    meanlog10Eta = mean(log10(vEta(vEta~=0)),'omitnan');
%     thresh_range = round(meanlog10Eta,0)-3:1:round(meanlog10Eta,0)+3;
    thresh_range = -5:0;
    xtextpos = log10T(60);
    hold on
    for i = thresh_range
        plot(log10T,i-log10T,'w--','LineWidth',1.5);
        text(xtextpos,i-xtextpos+0.4,num2str(i),'Color','w','FontWeight','b');
        xtextpos = xtextpos + 0.3;
    end
    hold off
%     ax1.XLim = [-6 -0.7];
%     ax1.YLim = [-2. 3.];
    title('(a)','FontSize',16)

    subplot(2,1,2)
    nbins = 40;
    plot_NND1Dhist(log10vEta,nbins)
    ax2 = gca;
%     ax2.YLim = [0 0.38];
%     ax2.XLim = [-9 0.5];
    hold on
    for i = thresh_range
        plot([i,i],[0,ax2.YLim(2)-0.03],...
                '--','Color',[0 0 0 0.5],'LineWidth',1.0);
        text(i,ax2.YLim(2)-0.02,num2str(i),'Color','k','FontWeight','b',...
            'HorizontalAlignment','center');
    end
    hold off
    
    title('(b)','FontSize',16)

end