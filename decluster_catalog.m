%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 14 October 2024
%   ...
%   version 1.2.0, 26 December 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.sRegion = 'italy';    % 'southcal'; % 
Model.fMc = 3.0;            % lower magnitude cutoff to perform analysis
sDeclstMethod = 'ML_NND';   % 'ETAS_SD';  % 'ZBZ2020';  % 'OriginalNND'; % 
sNNDmethod = 'OriginalNND'; % 'ZBZ2020';  % 
sDirPart = [Model.sRegion,'/'];
run([Model.sRegion,'/etas2d_parameters']);
nnd_parameters;
sDir_out = [Model.sRegion,'/decluster/'];
sFileName = [Model.sSeismName,'_',num2str(Model.fMc),'m_',num2str(Model.fTs),'d_',num2str(Model.fTe),'d_',sDeclstMethod];  % output file name
Model.sFileName = [sDir_out,sFileName];
% extracting earthquakes from a catalogue
tCat = readtable_catalog(Model.sEqCatName,'MagMin',Model.fMmin,'MagMax',Model.fMmax,'DateStart',Model.vDateStart,'Tstart',Model.fT0,'Tend',Model.fTe,...
                              'DepthMax',Model.fDepthMax,'SeismRegion',Model.vSeismReg_targ);
vCat = [tCat.Time, tCat.Lat, tCat.Lon, tCat.Mag, tCat.Depth]; % catalog in the target region
Model = set_etas_model(vCat,Model,PPip,'OptimMethod',sSolver);
%Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc)];
Model.sTitle = []; EqClassifier = []; tRealCatFeat = [];
if strcmp(sDeclstMethod,'OriginalNND')
    vCat = vCat(vCat(:,1) >= Model.fTs,:);    
elseif strcmp(sDeclstMethod,'ZBZ2020')
    %NNDpar.b = 0;
    NNDpar.M = 50;       % number of reshuflings of the background catalogue
    if strcmp(Model.sRegion,'southcal')
        NNDpar.alpha0 = 0.7; % \alpha_0 parameter
    elseif strcmp(Model.sRegion,'italy')
        NNDpar.alpha0 = 1.2; % \alpha_0 parameter
    end
    sNNDmethod = 'ZBZ2020';  % 
    vCat = vCat(vCat(:,1) >= Model.fTs,:);
elseif strcmp(sDeclstMethod,'ML_NND')
    sClassifierFile = [Model.sRegion,'/ml_nnd/',sFileName,'_classifier.mat'];   %  if already a trained classifier exists then load it from the file
    sFilePattern{1} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat.dat';
    sFilePattern{2} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat_par.dat';
    if isempty(sClassifierFile)
        sDir_train = [Model.sRegion,'/etas_train_cat'];
        [targetFiles, parFiles] = get_feature_file_names(sDir_train,sFilePattern);
        [tTrainCats, tTrainPars] = read_feature_files(targetFiles,parFiles,Model.fTs);
        EqClassifier = ml_nnd_training(tTrainCats,tTrainPars);
        %save(sClassifierFile_0,"EqClassifier");
    else
        load(sClassifierFile);
    end
    % read in the real catalog features
    sDir_fit = [Model.sRegion,'/etas_fit'];
    [targetFiles, parFiles] = get_feature_file_names(sDir_fit,sFilePattern);
    [tRealCatFeat,tRealPar] = read_feature_files(targetFiles,parFiles,Model.fTs);
    tRealCatFeat(any(isinf(tRealCatFeat{:,5:35}),2),:) = [];
    tRealCatFeat = inregion_events(tRealCatFeat,Model.vSeismReg_targ);
    vCat = [tRealCatFeat.Time, tRealCatFeat.Lat, tRealCatFeat.Lon, tRealCatFeat.Mag]; % catalog in the target region
elseif strcmp(sDeclstMethod,'ETAS_SD')
    % read in the real catalog features
    sDir_fit = [Model.sRegion,'/etas_fit'];
    sFilePattern{1} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat.dat';
    sFilePattern{2} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat_par.dat';
    [targetFiles, parFiles] = get_feature_file_names(sDir_fit,sFilePattern);
    [tRealCatFeat,tRealPar] = read_feature_files(targetFiles,parFiles,Model.fTs);
    tRealCatFeat(any(isinf(tRealCatFeat{:,5:35}),2),:) = [];
    tRealCatFeat = inregion_events(tRealCatFeat,Model.vSeismReg_targ);
    vCat = [tRealCatFeat.Time, tRealCatFeat.Lat, tRealCatFeat.Lon, tRealCatFeat.Mag]; % catalog in the target region
end
% declustering
[BkgrdInds, AshkInds] = event_decluster(vCat,NNDpar,Model,'NNDmethod',sNNDmethod,'DeclstMethod',sDeclstMethod,'Classifier',EqClassifier,'tCat',tRealCatFeat,'SaveFigure',Model.sFileName);
%
% frequency-magnitude statistics for earthquakes used to fit the PointProc model
Model.sFileName = [Model.sFileName,'_background'];
model_freq_mag(vCat(BkgrdInds,[1,4]),Model,Model.sTitle);
Model.sFileName = replace(Model.sFileName,'_background','_aftershocks');
model_freq_mag(vCat(AshkInds,[1,4]),Model,Model.sTitle);
