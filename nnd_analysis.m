%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 2.3.2, 2 December 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.sRegion = 'italy';    % 'southcal'; % 
Model.fMc = 3.0;            % lower magnitude cutoff to perform analysis
sNNDmethod = 'OriginalNND'; % 'ZBZ2020';     % 
sCatalog = 'real';          % 'synthetic';     % 
sDirPart = [Model.sRegion,'/'];
run([Model.sRegion,'/etas2d_parameters']);
nnd_parameters;
sDir_out = [Model.sRegion,'/nnd/'];
if strcmp(sCatalog,'synthetic')
    sDir_out = [Model.sRegion,'/nnd_synth/'];
    Model.sEqCatName  = [Model.sRegion,'/etas_test_cat/synth_southcal_catalog_150_feat.dat'];
    Model.sSeismName = 'synth_cat_150';
end
sFileName = [Model.sSeismName,'_',num2str(Model.fMc),'m_',num2str(Model.fTs),'d_',num2str(Model.fTe),'d_',sPointProc];  % output file name
Model.sFileName = [sDir_out,sFileName];
% extracting earthquakes from a catalogue
if strcmp(sCatalog,'real')
    tCat = readtable_catalog(Model.sEqCatName,'MagMin',Model.fMmin,'MagMax',Model.fMmax,'DateStart',Model.vDateStart,'Tstart',Model.fT0,'Tend',Model.fT1,...
                             'DepthMax',Model.fDepthMax,'SeismRegion',Model.vSeismReg_targ);
    vCat = [tCat.Time, tCat.Lat, tCat.Lon, tCat.Mag, tCat.Depth];
elseif strcmp(sCatalog,'synthetic')
    tCat = readtable(Model.sEqCatName);      % loading a user specified catalog file 
    % indxm = (tCat.Mag >= Model.fMmin) & (tCat.Mag <= Model.fMmax);
    % tCat = tCat(indxm,:);
    vCat = [tCat.Time, tCat.Lat, tCat.Lon, tCat.Mag];
end
vCat = vCat(vCat(:,1) >= Model.fTs,:);    
vTM = [tCat.Time, tCat.Mag];
if strcmp(sNNDmethod,'ZBZ2020')
    %NNDpar.b = 0;
    NNDpar.M = 50;       % number of reshuflings of the background catalogue
    if strcmp(Model.sRegion,'southcal')
        NNDpar.alpha0 = -1.0; % \alpha_0 parameter
    elseif strcmp(Model.sRegion,'italy')
        NNDpar.alpha0 = 1.0; % \alpha_0 parameter
    end
end
Model = set_etas_model(vCat,Model,PPip,'OptimMethod',sSolver);
%Model.sTitle = [Model.sTitle0,': T_s = ', num2str(Model.fTs), ', T_e = ', num2str(Model.fTe),'; N = ', num2str(Model.nNum), '; m_c = ',num2str(Model.fMc)];
Model.sTitle = [];
% nnd declustering
NND = model_nnd(vCat,NNDpar,Model,'NNDmethod',sNNDmethod);
%
% frequency-magnitude statistics for earthquakes used to fit the PointProc model
model_freq_mag(vTM(Model.nJs:Model.nJe,:),Model,Model.sTitle);
