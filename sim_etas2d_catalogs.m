%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 April 2021
%   ...
%   version: 2.1.0, 30 October 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.sRegion = 'italy';     % 'southcal'; % 
nTrainCats  = 100;           % number of training catalogues to simulate
nTestCats   = 100;           % number of testing catalogues to simulate
nCatalogs   = nTrainCats + nTestCats; % number of catalogues to simulate
sDirPart    = [Model.sRegion,'/'];
sPointProc  = 'ETAS2D8P';    % point process to use
sBackground = 'LoadFile';    %'Fractal'; %'Uniform'; %'LoadFile'; %'Catalog';  % background rate can be simulated, estimated from a catalogue or loaded from a file
Model.fMc   = 3.0;
run([Model.sRegion,'/etas2d_parameters']);
sDir_fit    = [Model.sRegion,'/etas_fit/'];
sDir_train  = [Model.sRegion,'/etas_train_cat/'];
sDir_test   = [Model.sRegion,'/etas_test_cat/'];
sFileName   = sprintf('synth_%s_catalog',Model.sRegion);
Model.sFileName = [sDir_train,sFileName];
% ETAS Parameters
vPar_est    = readmatrix([sDir_fit,Model.sRegion,'_catalog_feat_par.dat']);
Model.fBeta = vPar_est(1,9)*log(10);
if strcmp(sPointProc,'ETAS2D7P')
    sModPar = {'\mu', 'A',   '\alpha', 'c',   'p',   'd',   'q'};
    vPar    = vPar_est(1,1:7);           % for the model with 7 parameters
elseif strcmp(sPointProc,'ETAS2D8P')
    sModPar = {'\mu', 'A',   '\alpha', 'c',   'p',   'd',   'q',  '\gamma'};
    vPar    = vPar_est(1,1:8);
end
disp(['Branching ratio: ',num2str(vPar(2)*Model.fBeta/(Model.fBeta-vPar(3)))])        % n = A*beta/(beta-alpha)
global vEqCat;
global nEqNum;
nNmax = 100000;                     % the maximum number of events to generate
% load the background seismicity
if strcmp(sBackground,'LoadFile') % background is loaded from a file
    BackgrndRate = get_bckgrnd_rate(sBackground,[sDir_fit,Model.sRegion,'_catalog_est_background_rate.dat'],'Faults',[]);  % load the background rate u(x,y) from a file
    Model.vX = BackgrndRate.vX; % longitudes
    Model.vY = BackgrndRate.vY; % latitudes
    if strcmp(Model.sMapUnit,'degree')
        % convert to km because etas2d_simulator() generates the locations in km
        vXY_tr = coord_projection([Model.vY, Model.vX],'MapProjection',Model.sMapProj); % it returns [x, y] becasue the background rate is saved using longs and lats
        Model.vX = vXY_tr(:,1);
        Model.vY = vXY_tr(:,2);
    end
end
Model.mU  = BackgrndRate.mU;
Model.mMu = vPar(1)*BackgrndRate.mU;
nnd_parameters;
% simulate the model to generate a synthetic catalogue of earthquakes
for i = 1:(nTrainCats + nTestCats)
    vEqCat = zeros(nNmax,5);             % the earthquake catalogue
    nEqNum = 0;                          % the number of earthquakes generated: global variable
    if i <= nTrainCats
        disp(['Training catalog #',num2str(i)])
        Model.sFileName = [sDir_train,sFileName,'_',num2str(i)];
    else
        disp(['Testing catalog #',num2str(i-nTrainCats)])
        Model.sFileName = [sDir_test,sFileName,'_',num2str(i)];
    end
    etas2d_simulator(Model.fBeta,vPar,Model.fMmin,Model.fM0,Model.fT0,Model.fTe,nNmax,Model.vReg_targ,'PointProcess',sPointProc,'BinMag',Model.fDm,...
             'SaveCat',Model.sFileName,'Background',Model.vX,Model.vY,BackgrndRate.mU,'Anisotropy',bAnisotropy);
    vBgEvents_true = vEqCat(vEqCat(:,5) == 1,:);           % background events from the catalogue. Not applicable for real seismicity
    Model.sTitle = {[replace(Model.sSeismName,"_"," "),': $T_0 = ',num2str(Model.fT0),'$,  $T_e = ',num2str(Model.fTe),...
                '$;  $N_\mathrm{eq} = ',num2str(length(vEqCat(:,1))),'$; $m_c = ',num2str(Model.fMc),'$; $m_0 = ',num2str(Model.fM0),'$'],['$(',num2str(vPar,'%.4g '),')$']};
    Bckgrnd_est = etas2d_kernel_est_bckgrnd(Model.fTe,vEqCat,Model,vPar,'PointProcess',sPointProc);
    inR = inpolygon(vEqCat(:,3),vEqCat(:,2),Model.vReg_targ(:,1),Model.vReg_targ(:,2)); % logical indices 
    if strcmp(Model.sMapUnit,'degree')
        Model_degree = Model;
        vEqCat(:,2:3) = coord_projection(vEqCat(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % vCat = [lat, lon] input [x, y]
        vBgEvents_true(:,2:3) = coord_projection(vBgEvents_true(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % Bckgrnd_est.vBgDeclstr = [lat, lon] input [x, y]
        latlon = coord_projection([Model.vX, Model.vY(1)*ones(Model.nX,1)],'MapProjection',Model.sMapProj,'Direction','inverse'); % x-coordinate is changing with fixed y-coordinate
        Model_degree.vX = latlon(:,2); % longitude
        latlon = coord_projection([Model.vX(1)*ones(Model.nY,1), Model.vY],'MapProjection',Model.sMapProj,'Direction','inverse');  % y-coordinate is changing with fixed x-coordinate
        Model_degree.vY = latlon(:,1); % latitude

        vCat = vEqCat(inR,:);
        vU = Bckgrnd_est.vU(inR,:);
        vBgProb = Bckgrnd_est.vBgProb(inR,:);
        write_feature_catalog(Model.sFileName,Model,vCat,vPar,'Bkgrd_Int',vU,'Bkgrd_Prob',vBgProb);
        %write_feature_catalog(Model.sFileName,Model,vEqCat,vPar,'Bkgrd_Int',Bckgrnd_est.vU,'Bkgrd_Prob',Bckgrnd_est.vBgProb);
        cGeoFeatures.TargetRegion = coord_projection(Model.vReg_targ,'MapProjection',Model.sMapProj,'Direction','inverse');
        plot_etas2d_seismicity(vEqCat,Model_degree,'SaveFigure','TrueBckgrEvents',vBgEvents_true,'GeoFeatures',cGeoFeatures);
        close all
    end
end
