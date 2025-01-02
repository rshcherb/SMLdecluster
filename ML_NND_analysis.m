%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 9 August 2024
%   ...
%   version 1.8.0, 3 December 2024
%
clear all
addpath(genpath('MATLAB/'));
%%
Model.sRegion = 'italy';    % 'southcal'; % 
sNNDmethod    = 'ZBZ2020';  % 'OriginalNND'; % 
sClassifier   = 'Train';    % 'Load';     %   load or train classifier
Model.fMc = 3.0;            % lower magnitude cutoff to perform analysis
%
sDirPart = [Model.sRegion,'/'];
run([Model.sRegion,'/etas2d_parameters']);
nnd_parameters;
sDir_out = [Model.sRegion,'/ml_nnd/'];
Model.sFileName = [sDir_out,Model.sSeismName,'_',num2str(Model.fMc),'m_',num2str(Model.fTs),'d_',num2str(Model.fTe),'d_',sNNDmethod];  % output file name
sFilePattern{1} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat.dat';
sFilePattern{2} = 'catalog_' + optionalPattern(digitsPattern()+'_') + 'feat_par.dat';
sClassifierFile = [replace(Model.sFileName,sNNDmethod,''),'ML_NND_classifier.mat']; % file name for the classifier to save or load
%%
if strcmp(sClassifier,'Train')
    sDir_train = [Model.sRegion,'/etas_train_cat'];
    [targetFiles, parFiles] = get_feature_file_names(sDir_train,sFilePattern);
    [tTrainCats, tTrainPars] = read_feature_files(targetFiles,parFiles,Model.fTs);
    EqClassifier = ml_nnd_training(tTrainCats,tTrainPars);
    save(sClassifierFile,"EqClassifier");
elseif strcmp(sClassifier,'Load')
    load(sClassifierFile);
end
%%
fprintf('Testing catalog\n')
% create a list of the unique set of synthetic testing catalogs to read in
sDir_test = [Model.sRegion,'/etas_test_cat'];
[targetFiles, parFiles] = get_feature_file_names(sDir_test,sFilePattern);
[tTestCats, tTestPars] = read_feature_files(targetFiles,parFiles,Model.fTs);
% read in the real catalog
sDir_fit = [Model.sRegion,'/etas_fit'];
[targetFiles, parFiles] = get_feature_file_names(sDir_fit,sFilePattern);
[tRealCat, tRealPar] = read_feature_files(targetFiles,parFiles,Model.fTs);
tRealCat(any(isinf(tRealCat{:,5:35}),2),:) = [];

tTestCats = inregion_events(tTestCats,Model.vSeismReg_targ);
tRealCat = inregion_events(tRealCat,Model.vSeismReg_targ);

if strcmp(sNNDmethod,'ZBZ2020')
    %NNDpar.b = 0;
    NNDpar.M = 50;       % number of reshuflings of the background catalogue
    if strcmp(Model.sRegion,'southcal')
        NNDpar.alpha0 = 0.7; % \alpha_0 parameter
    elseif strcmp(Model.sRegion,'italy')
        NNDpar.alpha0 = 1.2; % \alpha_0 parameter
    end
end

%% NND analysis of the real and synthetic catalogs
vRealCat = [tRealCat.Time, tRealCat.Lat, tRealCat.Lon, tRealCat.Mag];
Model = set_etas_model(vRealCat,Model,PPip,'OptimMethod',sSolver);
[NND, NNDpar.gmmRealCat] = model_nnd(vRealCat,NNDpar,Model,'Display','off');
NNDpar.vLog10EtaThresh = NND.vLog10EtaThresh;
NNDpar.gmmTestCats = gmm_best(tTestCats.eta1,NNDpar.nGMMmax,'ThresholdType','minsaddle');

% Plot the NND 1D & 2D distributions
if strcmp(Model.sRegion,'southcal')
    NNDpar.EtaYLim = [0, 0.2];
    cFigLabel = {'(a)', '(c)'};
elseif strcmp(Model.sRegion,'italy')
    NNDpar.EtaYLim = [0, 0.3];
    cFigLabel = {'(b)', '(d)'};
end
plot_NND_T_R_Eta(tTestCats,tRealCat,NNDpar,'Region',Model.sRegion,'FigLabel',cFigLabel,'SaveFigure',Model.sFileName);
%% Plot the declustering comparison for the real catalog
if strcmp(Model.sRegion,'southcal')
    NNDpar.EtaYLim = [0, 0.75];
elseif strcmp(Model.sRegion,'italy')
    NNDpar.EtaYLim = [0, 0.65];
end
plot_NND_ML_real(tTestCats,tRealCat,NNDpar,Model,EqClassifier,'NNDmethod',sNNDmethod);
%% Plot the declustering comparison for the synthetic catalogs
plot_NND_ML_synth(tTestCats,tTestPars,NNDpar,Model,EqClassifier,'NNDmethod',sNNDmethod);

