function [BkgrdInds, AshkInds] = event_decluster(vCat,NNDpar,Model,options)
%
%   Perfomrs declustering of the catalog of events
%   vCat - catalogue of events
%   Model - structure that defines the seismicity model
%   OprArgs.FitModel - which mixture models to fit: 'GMM', 'Weilbul' 
%   OptArgs.FMD - FMD to use: 'GR', 'Exp', 'TapTapPareto', 'TapParetoPareto'
%   OptArgs.ThresholdType - 'intersect' the thresholds are computed as intersection of GMM compenets;
%                           'minsaddle' the thresholds are computed as the local minima for each saddle if present
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 14 October 2024
%   ...
%   version 1.5.0, 30 October 2024
%
    arguments
        vCat double
        NNDpar struct
        Model struct
        options.FMD char = 'GR'                  % frequency-magnitude distribution to rescale \eta, T, and R: 'GR', 
        options.DeclstMethod char = 'OriginalNND'  % which NND declustering method to use: 'OriginalNND', 'ZBZ2020', 'ML_NND', 'ETAS_SD'
        options.NNDmethod char = 'OriginalNND'   % which NND method to use: 'OriginalNND', 'ZBZ2020'
        options.Classifier struct = []
        options.tCat table = []
        options.TimeUnit char = 'day'
        options.Display char = 'on'              % 'on'/'off' to plot the distributions or not 
        options.SaveFigure char = []
        %OptArgs.BridgeTest char = 'on'           % 'on'/'off' to perform the Brigde test or not 
    end
    vPos = [400,300,800,600];
    vPos2pan = [400,300,800,800];
    
    % all declustering methods except ZBZ2020 use the original NND method when performing NND analysis
    [NND, ~, ~, ~, ~] = model_nnd(vCat,NNDpar,Model,'NNDmethod',options.NNDmethod,'Display','off');
    NNDpar.vLog10EtaThresh = NND.vLog10EtaThresh;

    if strcmp(options.DeclstMethod,'OriginalNND')
        BkgrdInds = NND.vClust == 0; % identify background events
        AshkInds = ~BkgrdInds;
        crossThreshBkgrdInds = [];
        crossThreshAshkInds = [];
    elseif strcmp(options.DeclstMethod,'ZBZ2020')
        BkgrdInds = NND.vClust == 0; % identify background events
        AshkInds = ~BkgrdInds;
        crossThreshBkgrdInds = BkgrdInds & (NND.vEta <= NND.vEtaThresh);
        crossThreshAshkInds = AshkInds & (NND.vEta > NND.vEtaThresh);
        disp(['alpha0 = ',num2str(NNDpar.alpha0)])
    elseif strcmp(options.DeclstMethod,'ML_NND')
        %Classifier = ml_nnd_declustering(tTrainCats,tTrainPars,Model.sRegion);
        eqClass = options.Classifier.predictFcn(options.tCat);
        AshkInds = eqClass == 1;
        BkgrdInds = eqClass == 0;
        crossThreshBkgrdInds = BkgrdInds & (NND.vEta <= NND.vEtaThresh);
        crossThreshAshkInds = AshkInds & (NND.vEta > NND.vEtaThresh);
        vCat = [options.tCat.Time, options.tCat.Lat, options.tCat.Lon, options.tCat.Mag];
    elseif strcmp(options.DeclstMethod,'ETAS_SD')
        U = rand(size(options.tCat.Bkgrd_Prob));
        BkgrdInds = U < options.tCat.Bkgrd_Prob;
        AshkInds = ~BkgrdInds;
        crossThreshBkgrdInds = BkgrdInds & (NND.vEta <= NND.vEtaThresh);
        crossThreshAshkInds = AshkInds & (NND.vEta > NND.vEtaThresh);
    end

    % perform the statistical tests
    vTimes = vCat(BkgrdInds,1);
    %vTimes = vCat(AshkInds,1);
    if 0    % synthetic Poisson times
        tmin = min(vCat(BkgrdInds,1));
        tmax = max(vCat(BkgrdInds,1));
        vTimes = sort(tmin + (tmax-tmin)*rand(length(vCat(BkgrdInds,1)),1));  % generate N random times in the interval [0, T]
        vCat(BkgrdInds,1) = vTimes;
    end

    % tic
    % [TestResults.Bridge.h, TestResults.Bridge.p, TestResults.Bridge.tstat] = bridge_test(vTimes,'Display',options.Display);
    % toc

    tic
    [TestResults.KS.h, TestResults.KS.p, TestResults.KS.tstat] = ks_test(vTimes,'Display',options.Display);
    toc

    tic
    [TestResults.BZ.h, TestResults.BZ.p, TestResults.BZ.tstat] = bz_test(vTimes,'Display',options.Display);
    toc
    
    if strcmp(options.Display,'on')
        figure('Name','T and R: background and aftershocks','Position',vPos);
        plot_TR_declustered(log10(NND.vT),log10(NND.vR),BkgrdInds,AshkInds,crossThreshBkgrdInds,crossThreshAshkInds,NNDpar);
        if strcmp(options.DeclstMethod,'ZBZ2020')
            text(0.05,0.1,['$\alpha_0 = ',num2str(NNDpar.alpha0),'$'],'Units','normalized','Interpreter','latex','FontSize',12);
        end
        if ~isempty(options.SaveFigure)
            save_cf(gcf,[options.SaveFigure,'_TR_declustered'],'fig','png','pdf');
        end

        figure('Name','Eta: background and aftershocks','Position',vPos);
        plot_eta_declustered(log10(NND.vEta),BkgrdInds,AshkInds,NNDpar,'PDF',NNDpar.sEtaDistrbFit,'PieSize',0.15);
        if ~isempty(options.SaveFigure)
            save_cf(gcf,[options.SaveFigure,'_Eta_declustered'],'fig','png','pdf');
        end

        figure('Name','Cumulative rate: background and aftershocks','Position',vPos);
        plot_cumul_rates(vCat,BkgrdInds,AshkInds,Model.vDateStart,'PieSize',0.15,'TestResults',TestResults);
        if ~isempty(options.SaveFigure)
            save_cf(gcf,[options.SaveFigure,'_cumulative_rates'],'fig','png','pdf');
        end

        figure('Name','Events: background and aftershocks','Position',vPos2pan);
        plot_events_time(vCat,BkgrdInds,AshkInds,Model);
        if ~isempty(options.SaveFigure)
            save_cf(gcf,[options.SaveFigure,'_events_time'],'fig','png','pdf');
        end
    end
end

