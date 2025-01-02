function plot_NND_ML_real(tTestCats,tRealCat,NNDpar,Model,EqClassifier,options)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 21 October 2024
%   ...
%   version 1.4.1, 6 December 2024
%
    arguments
        tTestCats table
        tRealCat table
        NNDpar struct
        Model struct
        EqClassifier struct
        options.SaveFig logical = true
        options.NNDmethod char = 'OriginalNND'   % which NND method to use: 'OriginalNND', 'ZBZ2020'
        options.FigLabel cell = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'} % letter lables for each figure
    end

    % identify NND model predictions based on the \eta thresholds
    if strcmp(options.NNDmethod,'OriginalNND')
        DclstNND = model_classification(tRealCat,NNDpar,Model,'DeclstMethod','OriginalNND','Catalog','real');
    elseif strcmp(options.NNDmethod,'ZBZ2020')
        DclstNND = model_classification(tRealCat,NNDpar,Model,'DeclstMethod','ZBZ2020','Catalog','real');
    end
    
    % identify trainedModel predictions. Use EqClassifier to make predictions on feature catalogs
    %DclstMLNNDsynth = model_classification(tTestCats,NNDpar,Model,'DeclstMethod','ML_NND','Classifier',EqClassifier,'Catalog','synthetic');
    DclstMLNNDreal  = model_classification(tRealCat,NNDpar,Model,'DeclstMethod','ML_NND','Classifier',EqClassifier,'Catalog','real');
    
    % identify ETAS SD model predictions
    %DclstETASSDsynth = model_classification(tTestCats,NNDpar,Model,'DeclstMethod','ETAS_SD','Catalog','synthetic');
    DclstETASSDreal  = model_classification(tRealCat,NNDpar,Model,'DeclstMethod','ETAS_SD','Catalog','real');
    
    trainingFeatures = print_features(EqClassifier.Features);
    fprintf('%s trained on %d synth catalogs with Features: %s\n\n',...
        EqClassifier.type,length(EqClassifier.trainingCatalogs),erase(trainingFeatures,"\"))

    % set plotting visuals
    TitleFontSize = 14;
    piesize = 0.11;
    xlbl = -0.12; ylbl = 1.05;
    
    % plot results
    f1 = figure('Name',[EqClassifier.type,' Real Catalog Comparison'],'Position',[400 300 900 600]);
    t = tiledlayout(f1,8,29,'TileSpacing','compact','Padding','compact');

    % plot NND distribution
    TRAx(1) = nexttile(t,1,[4,9]);
    if strcmp(options.NNDmethod,'OriginalNND')
        plot_TR_declustered(tRealCat.T1,tRealCat.R1,DclstNND.BkgrdInds,DclstNND.AshkInds,DclstNND.crossThreshBkgrdInds,...
            DclstNND.crossThreshAshkInds,NNDpar,'Legend','off','ThresholdValue','off');
        title('Original NND','FontSize',TitleFontSize,'Interpreter','latex');
    elseif strcmp(options.NNDmethod,'ZBZ2020')
        plot_TR_declustered(DclstNND.log10vT,DclstNND.log10vR,DclstNND.BkgrdInds,DclstNND.AshkInds,DclstNND.crossThreshBkgrdInds,...
            DclstNND.crossThreshAshkInds,NNDpar,'Legend','off','ThresholdValue','off');
        title('ZB2020 NND','FontSize',TitleFontSize,'Interpreter','latex');
    end
    text(TRAx(1),xlbl,ylbl,options.FigLabel{1},'FontSize',15,'Units','normalized','Interpreter','latex');
    % plot trained Model distribution
    TRAx(2) = nexttile(t,11,[4,9]);
    plot_TR_declustered(tRealCat.T1,tRealCat.R1,DclstMLNNDreal.BkgrdInds,DclstMLNNDreal.AshkInds,DclstMLNNDreal.crossThreshBkgrdInds,...
        DclstMLNNDreal.crossThreshAshkInds,NNDpar,'Legend','off','YAxisLabel','off','ThresholdValue','off');
    title('SML Method','FontSize',TitleFontSize,'Interpreter','latex');
    text(TRAx(2),xlbl,ylbl,options.FigLabel{2},'FontSize',15,'Units','normalized','Interpreter','latex');
    % plot ETAS SD Model distribution
    TRAx(3) = nexttile(t,21,[4,9]);
    plot_TR_declustered(tRealCat.T1,tRealCat.R1,DclstETASSDreal.BkgrdInds,DclstETASSDreal.AshkInds,DclstETASSDreal.crossThreshBkgrdInds,...
        DclstETASSDreal.crossThreshAshkInds,NNDpar,'Legend','off','YAxisLabel','off','ThresholdValue','off');
    title('Stochastic Declustering','FontSize',TitleFontSize,'Interpreter','latex');
    text(TRAx(3),xlbl,ylbl,options.FigLabel{3},'FontSize',15,'Units','normalized','Interpreter','latex');
    linkaxes(TRAx,'xy');

    % middle set of histograms for \eta
    hAx(1) = nexttile(t,117,[4,9]); % when using 8x29 grids
    %hAx(1) = nexttile(t,146,[4,9]);
    if strcmp(options.NNDmethod,'OriginalNND')
        plot_eta_declustered(tRealCat.eta1,DclstNND.BkgrdInds,DclstNND.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','PieSize',piesize);
    elseif strcmp(options.NNDmethod,'ZBZ2020')
        plot_eta_declustered(DclstNND.log10vEta,DclstNND.BkgrdInds,DclstNND.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','PieSize',piesize);
    end
    title(hAx(1),' ','FontSize',15); % this adds extra space between panels
    text(hAx(1),xlbl,ylbl,options.FigLabel{4},'FontSize',15,'Units','normalized','Interpreter','latex');
    hAx(2) = nexttile(t,127,[4,9]);
    plot_eta_declustered(tRealCat.eta1,DclstMLNNDreal.BkgrdInds,DclstMLNNDreal.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','YAxisLabel','off','PieSize',piesize);
    text(hAx(2),xlbl,ylbl,options.FigLabel{5},'FontSize',15,'Units','normalized','Interpreter','latex');
    title(hAx(2),' ','FontSize',15); % this adds extra space between panels
    hAx(3) = nexttile(t,137,[4,9]);
    plot_eta_declustered(tRealCat.eta1,DclstETASSDreal.BkgrdInds,DclstETASSDreal.AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'Legend','off','YAxisLabel','off','PieSize',piesize);
    text(hAx(3),xlbl,ylbl,options.FigLabel{6},'FontSize',15,'Units','normalized','Interpreter','latex');
    title(hAx(3),' ','FontSize',15); % this adds extra space between panels
    linkaxes(hAx,'xy');

    if options.SaveFig
        save_cf(gcf,[Model.sFileName,'_NND_ML_SD_real_catalog'],'fig','png');
    end
end