function plot_NND_ML_synth(tTestCats,ETASpar,NNDpar,Model,Classifier,options)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 October 2024
%   ...
%   version 1.4.0, 3 December 2024
%
    arguments
        tTestCats table
        ETASpar table
        NNDpar struct
        Model struct
        Classifier struct
        options.N_NND (1,1) = 1
        options.plotCat logical = true
        options.SaveFig logical = true
        options.NNDmethod char = 'OriginalNND'   % which NND method to use: 'OriginalNND' or 'ZBZ2020'
        options.FigLabel cell = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'} % letter lables for each figure
    end
    
    % identify NND model predictions
    DclstNND = model_classification(tTestCats,NNDpar,Model,'DeclstMethod',options.NNDmethod,'Catalog','synthetic');
    % identify Classifier predictions
    DclstMLNND = model_classification(tTestCats,NNDpar,Model,'DeclstMethod','ML_NND','Classifier',Classifier,'Catalog','synthetic');
    % identify SD_ETAS model predictions
    DclstETASSD = model_classification(tTestCats,NNDpar,Model,'DeclstMethod','ETAS_SD','Catalog','synthetic');

    % Set plotting visuals
    clrs = lines(5);
    colors = [clrs(1,:);clrs(2,:);clrs(4,:);clrs(3,:);[0 0 0]];
    txtSize = 11;
    LabelFontSize = 12;
    TitleFontSize = 14;
    xtext = 0.02; ytext = 0.12; dy = 0.07;
    xlbl = -0.14; ylbl = 1.05;

    % plot results
    f1 = figure('Name',[Classifier.type,' Synthetic Comparison'],'Position',[400 300 900 800]);
    t = tiledlayout(f1,13,15,'TileSpacing','Compact','Padding','Compact');

    % plot NND TR distribution
    TRAx(1) = nexttile(t,1,[5,5]);
    if strcmp(options.NNDmethod,'OriginalNND')
        plot_TR_declustered(tTestCats.("T"+num2str(options.N_NND)),tTestCats.("R"+num2str(options.N_NND)),...
                            DclstNND.BkgrdIndsWoutMissed,DclstNND.AshkIndsWoutMissed,DclstNND.MissedBkgrdInds,DclstNND.MissedAshkInds,NNDpar,...
                            'Legend','off','Colors',colors,'MarkerSize',5,'Synthetic',true,'ThresholdValue','off');
        title('Original NND','FontSize',TitleFontSize,'Interpreter','latex');
    elseif strcmp(options.NNDmethod,'ZBZ2020')
        plot_TR_declustered(DclstNND.log10vT,DclstNND.log10vR,...
                            DclstNND.BkgrdIndsWoutMissed,DclstNND.AshkIndsWoutMissed,DclstNND.MissedBkgrdInds,DclstNND.MissedAshkInds,NNDpar,...
                            'Legend','off','Colors',colors,'MarkerSize',5,'Synthetic',true,'ThresholdValue','off');
        title('ZB2020 NND','FontSize',TitleFontSize,'Interpreter','latex');
    end
    text(xtext,ytext-dy,"\textbf{Accuracy: "+num2str(DclstNND.Accuracy,'%.1f')+"\%}",'units','normalized','FontSize',txtSize,'HorizontalAlignment','l','Interpreter','Latex');
    text(TRAx(1),xlbl,ylbl,options.FigLabel{1},'FontSize',15,'Units','normalized','Interpreter','latex');
        
    % plot trained Model TR distribution
    TRAx(2) = nexttile(t,6,[5,5]);
    plot_TR_declustered(tTestCats.("T"+num2str(options.N_NND)),tTestCats.("R"+num2str(options.N_NND)),...
                        DclstMLNND.BkgrdIndsWoutMissed,DclstMLNND.AshkIndsWoutMissed,DclstMLNND.MissedBkgrdInds,DclstMLNND.MissedAshkInds,NNDpar,...
                        'Legend','off','Colors',colors,'MarkerSize',5,'Synthetic',true,'ThresholdValue','off');
    text(xtext,ytext-dy,"\textbf{Accuracy: "+num2str(DclstMLNND.Accuracy,'%.1f')+"\%}",'units','normalized','FontSize',txtSize,'HorizontalAlignment','l','Interpreter','Latex');
    text(TRAx(2),xlbl,ylbl,options.FigLabel{2},'FontSize',15,'Units','normalized','Interpreter','latex');
    title('SML Method','FontSize',TitleFontSize,'Interpreter','latex')
        
    % plot ETAS SD Model TR distribution
    TRAx(3) = nexttile(t,11,[5,5]);
    plot_TR_declustered(tTestCats.("T"+num2str(options.N_NND)),tTestCats.("R"+num2str(options.N_NND)),...
                        DclstETASSD.BkgrdIndsWoutMissed,DclstETASSD.AshkIndsWoutMissed,DclstETASSD.MissedBkgrdInds,DclstETASSD.MissedAshkInds,NNDpar,...
                        'Legend','off','Colors',colors,'MarkerSize',5,'Synthetic',true,'ThresholdValue','off');
    text(xtext,ytext-dy,"\textbf{Accuracy: "+num2str(DclstETASSD.Accuracy,'%.1f')+"\%}",'units','normalized','FontSize',txtSize,'HorizontalAlignment','l','Interpreter','Latex');
    text(TRAx(3),xlbl,ylbl,options.FigLabel{3},'FontSize',15,'Units','normalized','Interpreter','latex');
    title('Stochastic Declustering','FontSize',TitleFontSize,'Interpreter','latex')
    linkaxes(TRAx,'xy')
        
    % plot NND eta histogram
    hAx(1) = nexttile(t,76,[5,5]);
    if strcmp(options.NNDmethod,'OriginalNND')
        plot_eta_declustered(tTestCats.eta1,DclstNND.BkgrdInds,DclstNND.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','Title',' ');
    elseif strcmp(options.NNDmethod,'ZBZ2020')
        plot_eta_declustered(DclstNND.log10vEta,DclstNND.BkgrdInds,DclstNND.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','Title',' ');
    end
    text(hAx(1),xlbl,ylbl,options.FigLabel{4},'FontSize',15,'Units','normalized','Interpreter','latex');
    % plot trained Model eta histogram
    hAx(2) = nexttile(t,81,[5,5]);
    plot_eta_declustered(tTestCats.eta1,DclstMLNND.BkgrdInds,DclstMLNND.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','YAxisLabel','off','Title',' ');
    text(hAx(2),xlbl,ylbl,options.FigLabel{5},'FontSize',15,'Units','normalized','Interpreter','latex');
    % plot ETAS SD Model eta histogram
    hAx(3) = nexttile(t,86,[5,5]);
    plot_eta_declustered(tTestCats.eta1,DclstETASSD.BkgrdInds,DclstETASSD.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','YAxisLabel','off','Title',' ');
    text(hAx(3),xlbl,ylbl,options.FigLabel{6},'FontSize',15,'Units','normalized','Interpreter','latex');
    linkaxes(hAx,'xy')

    % plot NND confusion chart
    ccAx(1) = nexttile(t,152,[3,3]);
    plot_confusion_chart(tTestCats.Type2,DclstNND.Predict)
    xlabel('Predicted Class','FontSize',LabelFontSize,'Interpreter','latex');
    ylabel('True Class','FontSize',LabelFontSize,'Interpreter','latex');
    title(' ','FontSize',15); % this adds extra space between panels
    text(ccAx(1),xlbl,ylbl,options.FigLabel{7},'FontSize',15,'Units','normalized','Interpreter','latex');
        
    % plot trained Model confusion chart
    ccAx(2) = nexttile(t,157,[3,3]);
    plot_confusion_chart(tTestCats.Type2,DclstMLNND.Predict);
    xlabel('Predicted Class','FontSize',LabelFontSize,'Interpreter','latex');
    text(ccAx(2),xlbl,ylbl,options.FigLabel{8},'FontSize',15,'Units','normalized','Interpreter','latex');
        
    % plot ETAS SD Model confusion chart
    ccAx(3) = nexttile(t,162,[3,3]);
    plot_confusion_chart(tTestCats.Type2,DclstETASSD.Predict);
    xlabel('Predicted Class','FontSize',LabelFontSize,'Interpreter','latex');
    text(ccAx(3),xlbl,ylbl,options.FigLabel{9},'FontSize',15,'Units','normalized','Interpreter','latex');

    if options.SaveFig
        save_cf(gcf,[Model.sFileName,'_NND_ML_SD_synthetic_catalogs'],'fig','png');
    end    
end
