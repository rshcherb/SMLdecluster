%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 6 May 2021
%   ...
%   version 1.4.1, 19 November 2024
%
clear all
addpath(genpath('MATLAB/'));

Model.sRegion = 'italy';     % 'southcal';  % 
sPointProc  = 'ETAS2D8P';    % point process to use
Model.fMc   = 3.0;           % lower magnitude cutoff to perform analysis
sBackground = 'Catalog';     %'Fractal';     %'Uniform';  %  
sDirPart    = [Model.sRegion,'/'];
run([Model.sRegion,'/etas2d_parameters']);
sFileName = [Model.sRegion,'_catalog'];
Model.sFileName = [sDir_out,sFileName];
load_data; % extracting earthquakes from a catalogue
if strcmp(sBackground,'Catalog')
    Model.nX = 256;
    Model.nY = 256;
    Model.vX = linspace(min(Model.vReg_all(:,1)),max(Model.vReg_all(:,1)),Model.nX)';  % the x-coordinates of the background grid
    Model.vY = linspace(min(Model.vReg_all(:,2)),max(Model.vReg_all(:,2)),Model.nY)';  % the y-coordinates of the background grid
    %Bckgrnd0 = get_bckgrnd_rate('Uniform',Model,0.9e-6);  % uniform background rate fNormFactor*ones(x,y)
    Backgrnd0  = get_bckgrnd_rate('Catalog',vCat_km(vCat_km(:,1) <= Model.fTe,:),Model,'All');  % estimate the background rate from the catalogue
end
Model.mU = Backgrnd0.mU;  % set the initial background rate
% set the structure for the fit
Model = set_etas2d_model(vCat_km,Model,PPip,'Solver',sSolver);
Model.sTitle = [Model.sTitle0,'$T_s = ', num2str(Model.fTs),'$, $T_e = ',num2str(Model.fTe),'$; $N_\mathrm{targ} = ',num2str(Model.ETAS.nNtarg),'$; $m_c = ',num2str(Model.fMc),'$'];

[vPar_est, vParErr, fLle, Backgrnd_est] = model_etas2d_rate(vCat_km,Model,'PointProcess',sPointProc,'Background',sBackground);
plot_etas2d_rate(vCat_km,Model,vPar_est,vParErr,fLle,Backgrnd_est,'GeoFeatures',cGeoFeatures);

nnd_parameters;
vCat = vCat(Model.ETAS.inR,:);
vU = Backgrnd_est.vU(Model.ETAS.inR,:);
vBgProb = Backgrnd_est.vBgProb(Model.ETAS.inR,:);
write_feature_catalog(Model.sFileName,Model,vCat,vPar_est,vParErr,'N_NND',10,'Bkgrd_Int',vU,'Bkgrd_Prob',vBgProb);
% frequncy-magnitude statistics for earthquakes used to fit the PointProc model
% model_freq_mag(vCat(Model.ETAS.inTR,[1,4]),Model,Model.sTitle);
