function plot_NND_T_R_Eta(testCats,realCat,NNDpar,OptArgs)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 13 September 2024
%   ...
%   version 1.1.0, 3 December 2024
%
    arguments
        testCats table
        realCat table
        NNDpar struct
        OptArgs.SaveFigure char = []
        OptArgs.Region char = ''
        OptArgs.FigLabel cell = {'(a)', '(b)'} % letter lables for each figure
    end
    xlbl = -0.08; ylbl = 1.05;
    
    figure('Name',sprintf('NND Plots for %s',OptArgs.Region),'Position',[400 300 600 800]);
    t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
    
    nexttile(t);
    plot_NND_T_R(testCats,realCat,NNDpar)
    text(xlbl,ylbl,OptArgs.FigLabel{1},'FontSize',15,'Units','normalized','Interpreter','latex');
    
    nexttile(t);
    plot_NND_Eta(testCats,realCat,NNDpar)
    text(xlbl,ylbl,OptArgs.FigLabel{2},'FontSize',15,'Units','normalized','Interpreter','latex');
    
    if ~isempty(OptArgs.SaveFigure)
        save_cf(gcf,[OptArgs.SaveFigure,'_T_R_Eta_distributions'],'fig','png');
    end
end
