function Declust = model_classification(tFeatCat,NNDpar,Model,options)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 24 October 2024
%   ...
%   version 1.2.1, 31 October 2024
%
    arguments
        tFeatCat table
        NNDpar struct
        Model struct
        options.DeclstMethod char = 'OriginalNND'  % which NND declustering method to use: 'OriginalNND', 'ZBZ2020', 'ML_NND', 'ETAS_SD'
        options.Catalog char = 'real'              % type of data: 'real', 'synthetic'
        options.Classifier struct = []             % classfier structure for ML_NND method
        options.NNDFeat double = 0                 % whether to use NND features from the catalog tFeatCat or to compute
    end

    if strcmp(options.Catalog,'synthetic')
        TrueAshks = (tFeatCat.Type2 == 1);  % true aftershock indices from the synthetic catalog
        TrueBkgrd = (tFeatCat.Type2 == 0);  % true background event indices from the synthetic catalog
    end
    % Identify NND model predictions
    if strcmp(options.DeclstMethod,'OriginalNND')
        AshkInds = tFeatCat.eta1 <= NNDpar.vLog10EtaThresh;
        BkgrdInds = tFeatCat.eta1 > NNDpar.vLog10EtaThresh;
        Predict = AshkInds;
        crossThreshBkgrdInds = [];
        crossThreshAshkInds = [];
    elseif strcmp(options.DeclstMethod,'ZBZ2020')
        uniqueCats = unique(tFeatCat.Catalog); % indentify unique catalogue numbers in the tFeatCat
        nCat = length(uniqueCats);
        ntot = height(tFeatCat);
        BkgrdInds = false(ntot,1); % for background events
        AshkInds = false(ntot,1);  % for aftershocks
        crossThreshBkgrdInds = false(ntot,1);
        crossThreshAshkInds = false(ntot,1);
        if strcmp(options.Catalog,'synthetic')
            TrueAshks = false(ntot,1);  % true aftershock indices from the synthetic catalog
            TrueBkgrd = false(ntot,1);  % true background event indices from the synthetic catalog
        end
        Declust.log10vT = zeros(ntot,1);
        Declust.log10vR = zeros(ntot,1);
        Declust.log10vEta = zeros(ntot,1);
        for i = 1:nCat
            inds = (tFeatCat.Catalog == uniqueCats(i)); % extract only records for a given catalog number
            vCat = [tFeatCat.Time(inds), tFeatCat.Lat(inds), tFeatCat.Lon(inds), tFeatCat.Mag(inds)];            
            [NND, ~, ~, ~, ~] = model_nnd(vCat,NNDpar,Model,'NNDmethod','ZBZ2020','Display','off');
            Declust.log10vT(inds) = log10(NND.vT);
            Declust.log10vR(inds) = log10(NND.vR);
            Declust.log10vEta(inds) = log10(NND.vEta);
            %disp(['cat# ',num2str(uniqueCats(i)),', indx = ',num2str(sum(inds)),', vCat = ',num2str(size(vCat,1)),', vEta = ',num2str(size(NND.vEta,1))])
            disp(['cat# ',num2str(uniqueCats(i)),', N = ',num2str(sum(inds))])
            BkgrdInds(inds) = NND.vClust == 0; % identify background events
            AshkInds(inds) = ~BkgrdInds(inds);
            crossThreshBkgrdInds(inds) = BkgrdInds(inds) & (NND.vEta <= NND.vEtaThresh);
            crossThreshAshkInds(inds) = AshkInds(inds) & (NND.vEta > NND.vEtaThresh);
            if options.NNDFeat
                crossThreshBkgrdInds(inds) = BkgrdInds(inds) & (tFeatCat.eta1(inds) <= NND.vLog10EtaThresh);
                crossThreshAshkInds(inds) = AshkInds(inds) & (tFeatCat.eta1(inds) > NND.vLog10EtaThresh);
            end
            if strcmp(options.Catalog,'synthetic')
                TrueAshks(inds) = (tFeatCat.Type2(inds) == 1);  % true aftershock indices from the synthetic catalog
                TrueBkgrd(inds) = (tFeatCat.Type2(inds) == 0);  % true background event indices from the synthetic catalog
            end
        end
        Predict = AshkInds;
    elseif strcmp(options.DeclstMethod,'ML_NND')
        % Use Classifier to make predictions on tFeatCat
        Predict = options.Classifier.predictFcn(tFeatCat);
        AshkInds = (Predict == 1);
        BkgrdInds = (Predict == 0);
        crossThreshBkgrdInds = BkgrdInds & (tFeatCat.eta1 <= NNDpar.vLog10EtaThresh);
        crossThreshAshkInds = AshkInds & (tFeatCat.eta1 > NNDpar.vLog10EtaThresh);
    elseif strcmp(options.DeclstMethod,'ETAS_SD')
        % Identify Zhuang model predictions
        U = rand(size(tFeatCat.Bkgrd_Prob));
        BkgrdInds = U < tFeatCat.Bkgrd_Prob;
        AshkInds = ~BkgrdInds;
        Predict = AshkInds;
        crossThreshBkgrdInds = BkgrdInds & (tFeatCat.eta1 <= NNDpar.vLog10EtaThresh);
        crossThreshAshkInds = AshkInds & (tFeatCat.eta1 > NNDpar.vLog10EtaThresh);
    end

    % make a structure to return the values
    Declust.AshkInds  = AshkInds;   % aftershock indices: 1 if aftershock 0 otherwise
    Declust.BkgrdInds = BkgrdInds;  % background indices: 1 if background 0 otherwise
    Declust.Predict   = Predict;    % predicted indices: 0 for background and 1 for aftershocks
    if strcmp(options.Catalog,'real')
        Declust.crossThreshBkgrdInds = crossThreshBkgrdInds;
        Declust.crossThreshAshkInds  = crossThreshAshkInds;
    end
    if strcmp(options.Catalog,'synthetic')
        MissedAshkInds  = BkgrdInds & TrueAshks;
        MissedBkgrdInds = AshkInds & TrueBkgrd;
        CorrectInds = ~MissedAshkInds & ~MissedBkgrdInds;
        % Calculate model error and accuracy percentages
        Declust.AshkErr  = 100*(sum(MissedAshkInds)/sum(TrueBkgrd));
        Declust.BkgrdErr = 100*(sum(MissedBkgrdInds)/sum(TrueAshks));
        Declust.Accuracy = 100*(sum(CorrectInds)/size(tFeatCat,1));
        Declust.BkgrdIndsWoutMissed = TrueBkgrd & ~MissedBkgrdInds;
        Declust.AshkIndsWoutMissed  = TrueAshks & ~MissedAshkInds;
        % 
        Declust.MissedAshkInds  = MissedAshkInds;
        Declust.MissedBkgrdInds = MissedBkgrdInds;
        Declust.CorrectInds     = CorrectInds;
    end
end
