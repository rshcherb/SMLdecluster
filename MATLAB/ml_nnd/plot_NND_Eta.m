function plot_NND_Eta(testCats,realCat,NNDpar)
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
    txtSize = 10;
    EtahistXLim = [log10(NNDpar.EtaLim) 2];
    EtahistYLim = [NNDpar.EtaYLim 0.1];
    LabelFontSize = 14;
    LgdFontSize = 11;
    testCatsFaceAlpha = 0.5;
    realCatFaceAlpha = 0.6;

    % Set histogram properties
    nbins = 22;
    binWidth = 0.15;
    % Plot NND eta histogram
    ax = gca;
    colororder(ax,[colors(2,:);colors(1,:)])

    yyaxis left
    hTestCats = histogram(testCats.eta1,nbins,'Normalization','pdf',...
        'BinWidth',binWidth,'FaceColor',colors(2,:),'EdgeColor',colors(2,:),...
        'FaceAlpha',testCatsFaceAlpha,'DisplayName','Synthetic Events');
    ylabel('density','FontSize',LabelFontSize,'Interpreter','latex','Color','k');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.YLim = EtahistYLim(1:2);
    hold on

    yyaxis right
    hRealCat = histogram(realCat.eta1,nbins,'Normalization','pdf',...
        'BinWidth',binWidth*1.3,'FaceColor',colors(1,:),'EdgeColor',colors(1,:),...
        'FaceAlpha',realCatFaceAlpha,'DisplayName','Real Events');
    %ylabel('Real Event Count','FontSize',12,'Interpreter','latex');
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    plot([NNDpar.vLog10EtaThresh,NNDpar.vLog10EtaThresh],[0,EtahistYLim(2)],'--','Color',colors(3,:),'LineWidth',1.5,'DisplayName','NND Threshold');
    hold off
%     grid on
    
    [~,realCatClstIdx] = min(NNDpar.gmmRealCat.mu);
    [~,realCatBkgrdIdx] = max(NNDpar.gmmRealCat.mu);
    [~,testCatsClstIdx] = min(NNDpar.gmmTestCats.mu);
    [~,testCatsBkgrdIdx] = max(NNDpar.gmmTestCats.mu);
    ax.XLim = EtahistXLim(1:2);
    ax.YLim = EtahistYLim(1:2);
%     ax.XTick = round(ax.XLim(1),0):EtahistXLim(3):round(ax.XLim(2),0);
%     ax.YTick = round(ax.YLim(1),0):EtahistYLim(3):round(ax.YLim(2),0);
    %legend('Location','northeastoutside','box','off','FontSize',LgdFontSize,'Interpreter','latex')
    xlabel('Nearest-neighbour distance, $\log_{10} (\eta)$','FontSize',LabelFontSize,'Interpreter','latex')
    xtext = 0.13;
    ytext = 0.87;
    dy = 0.06;
    dx = 0.12;
    lblTextAlign = 'l';
    dataTextAlign = 'r';
    
    text(xtext + 0.05,ytext + 0.08,"\underline{\textbf{ GMM }}",...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',lblTextAlign,'Color',colors(3,:),'Interpreter','Latex');
    text(xtext,ytext + 0.01,"\textbf{Mean}",...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',lblTextAlign,'Color',colors(3,:),'Interpreter','Latex');
    text(xtext + dx,ytext + 0.01,"\textbf{Ratio}",...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',lblTextAlign,'Color',colors(3,:),'Interpreter','Latex');
    text(xtext - dx,ytext - 1.5*dy,"\textbf{Mode 1}",...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',lblTextAlign,'Color',colors(3,:),'Interpreter','Latex');
    text(xtext - dx,ytext - 4*dy,"\textbf{Mode 2}",...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',lblTextAlign,'Color',colors(3,:),'Interpreter','Latex');
    
    text(xtext + 0.5*dx,ytext - dy,sprintf('%.2f',NNDpar.gmmTestCats.mu(testCatsClstIdx)),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(2,:),'Interpreter','Latex');
    text(xtext + 0.5*dx,ytext - 2*dy,sprintf('%.2f',NNDpar.gmmRealCat.mu(realCatClstIdx)),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(1,:),'Interpreter','Latex');
    text(xtext + 0.5*dx,ytext - 3.5*dy,sprintf('%.2f',NNDpar.gmmTestCats.mu(testCatsBkgrdIdx)),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(2,:),'Interpreter','Latex');
    text(xtext + 0.5*dx,ytext - 4.5*dy,sprintf('%.2f',NNDpar.gmmRealCat.mu(realCatBkgrdIdx)),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(1,:),'Interpreter','Latex');
    
    text(xtext + 1.45*dx,ytext - dy,num2str(NNDpar.gmmTestCats.ComponentProportion(testCatsClstIdx),'%.2f'),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(2,:),'Interpreter','Latex');
    text(xtext + 1.45*dx,ytext - 2*dy,num2str(NNDpar.gmmRealCat.ComponentProportion(realCatClstIdx),'%.2f'),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(1,:),'Interpreter','Latex');
    text(xtext + 1.45*dx,ytext - 3.5*dy,num2str(NNDpar.gmmTestCats.ComponentProportion(testCatsBkgrdIdx),'%.2f'),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(2,:),'Interpreter','Latex');
    text(xtext + 1.45*dx,ytext - 4.5*dy,num2str(NNDpar.gmmRealCat.ComponentProportion(realCatBkgrdIdx),'%.2f'),...
        'units','normalized','FontSize',txtSize,...
        'HorizontalAlignment',dataTextAlign,'Color',colors(1,:),'Interpreter','Latex');
end