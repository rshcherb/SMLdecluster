function plot_etas2d_rate(vCat,Model,vPar_est,vParErr,fLle,Backgrnd_est,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 24 November 2023
%   ...
%   version 1.0.1, 26 November 2023
%
    sPointProc = 'ETAS2D8P';   % 
    bTrueBckgrndEvents = false;
    vBgEvents_true = [];
    cGeoFeatures = [];
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('TrueBckgrEvents',varargin{k})
            vBgEvents_true = varargin{k+1};
            if ~isempty(vBgEvents_true)
                bTrueBckgrndEvents = true;
            end
        end
        if strcmp('GeoFeatures',varargin{k})
            cGeoFeatures = varargin{k+1};
        end
    end
    vCat_m0 = vCat;
    vCat_m0(:,4) = vCat_m0(:,4) - Model.ETAS.fM0; % subtract the reference magnitude fM0
    
    [fMs, iMs] = max(vCat(:,4));      % the largest event in the catalogue
    fTplot     = vCat(iMs,1);         % time of the strongest event in days 
    fMsMin     = fMs - 1.0;           % strong events to plot

    Model.mMu = vPar_est(1)*Model.mU;
    Model.sTitle = {Model.sTitle,['(',num2str(vPar_est,'%.4g '),')']};

    % plot the ground intensity function at the time of the largest event
    %plot_etas2d_gif(fTplot,vCat,Model.vX,Model.vY,Bckgrnd_est.mMu,vPar_est,Model.fM0,Model.sTitle,'RateMax',Model.fRateMax,'MapUnits','km');
    plot_etas2d_gif(fTplot,vCat_m0,Model,Backgrnd_est.mMu,vPar_est,'SaveFigure',Model.sFileName,'GeoFeatures',cGeoFeatures);
    
    nXreg = Model.nX; %ceil(nX/16);
    nYreg = Model.nY; %ceil(nY/16);
    mSpatRes = etas2d_spat_residual(Model,Backgrnd_est.mMu,vCat_m0,vPar_est,Model.vReg_all,'DummyPoints','meshgrid',nXreg,nYreg);
    
    if strcmp(Model.sMapUnit,'degree')  % if true then convert into degrees
        vCat(:,2:3) = coord_projection(vCat(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % vCat = [lat, lon] input [x, y]
        %Backgrnd_est.vBgDeclstr(:,2:3) = coord_projection(Backgrnd_est.vBgDeclstr(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % Bckgrnd_est.vBgDeclstr = [lat, lon] input [x, y]
        if ~isempty(vBgEvents_true)
            vBgEvents_true(:,2:3) = coord_projection(vBgEvents_true(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % Bckgrnd_est.vBgDeclstr = [lat, lon] input [x, y]
        end
%         latlon = coord_projection(Model.vReg_all,'MapProjection',Model.sMapProj,'Direction','inverse');
%         Model.vReg_all = [latlon(:,2), latlon(:,1)]; % for plotting 
%         latlon = coord_projection(Model.vReg_targ,'MapProjection',Model.sMapProj,'Direction','inverse');
%         Model.vReg_targ = [latlon(:,2), latlon(:,1)]; % for plotting 
        latlon = coord_projection([Model.vX, Model.vY(1)*ones(Model.nX,1)],'MapProjection',Model.sMapProj,'Direction','inverse'); % x-coordinate is changing with fixed y-coordinate
        Model.vX = latlon(:,2); % longitude
        latlon = coord_projection([Model.vX(1)*ones(Model.nY,1), Model.vY],'MapProjection',Model.sMapProj,'Direction','inverse');  % y-coordinate is changing with fixed x-coordinate
        Model.vY = latlon(:,1); % latitude
    end
    % 
    plot_etas2d_seismicity(vCat,Model,'SaveFigure','TrueBckgrEvents',vBgEvents_true,'MsMagMin',fMsMin,'GeoFeatures',cGeoFeatures);
    if bTrueBckgrndEvents
        % plot the background rate, true background events and declustered events
        plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'TrueBckgrEvents',vBgEvents_true,'Declustered','GeoFeatures',cGeoFeatures);
    end
    % plot the estimated background rate
    plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'GeoFeatures',cGeoFeatures);
    % plot the background rate and declustered events
    plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'DeclusteredEvents','Dot','GeoFeatures',cGeoFeatures);
    % plot the background rate and declustered events
    plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'DeclusteredEvents','Circle','GeoFeatures',cGeoFeatures);
    % plot the background rate, declustered evenmts, and aftershocks
    plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'DeclusteredEvents','Dot','Aftershocks',vCat(~Backgrnd_est.indxBg,:),'GeoFeatures',cGeoFeatures);
    % plot the background rate, declustered evenmts, and aftershocks
    plot_etas2d_decluster(Backgrnd_est,Model,Backgrnd_est.mMu,'SaveFigure',Model.sFileName,'DeclusteredEvents','Circle','Aftershocks',vCat(~Backgrnd_est.indxBg,:),'GeoFeatures',cGeoFeatures);
    %
    % plot the results of the fit and stochastic declustering
    plot_etas2d_kernel_est_seism(vCat,Model,Backgrnd_est,'PointProcess',sPointProc,'SaveFigure',Model.sFileName,'TrueBckgrEvents',vBgEvents_true,'GeoFeatures',cGeoFeatures);
    % plot the diagnostics
    plot_etas2d_diagnostic(vCat,Model,vPar_est,vParErr,fLle,Backgrnd_est,mSpatRes,'SaveFigure',Model.sFileName,'GeoFeatures',cGeoFeatures);
end

