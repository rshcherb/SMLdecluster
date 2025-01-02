function plot_etas2d_diagnostic(vCat,Model,vETASPar,vParErr,fLle,Bckgrnd_est,mSpatRes,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 12 December 2021
%   ...
%   version 1.1.1, 18 October 2023
%
    sPointProc = 'ETAS2D8P';   % 
    bCoastline = false;
    bTargetReg = false;
    bSaveFig   = false;
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('GeoFeatures',varargin{k})
            cGeoFeatures = varargin{k+1};
            if isfield(cGeoFeatures,'Coastline')
                bCoastline = true;
            end
            if isfield(cGeoFeatures,'TargetRegion')
                bTargetReg = true;
            end
        end
        if strcmp('SaveFigure',varargin{k})
            bSaveFig = true;
            sFileName = varargin{k+1};
        end
    end
    vPos = [500 300 800 600];
    vPos2 = [500 300 700 600];
    % 
    fTs   = Model.ETAS.fTs;
    fTe   = Model.ETAS.fTe;
    if strcmp(sPointProc,'ETAS2D7P')
        vParName = {'\mu' 'A' '\alpha' 'c' 'p' 'd' 'q'};
    elseif strcmp(sPointProc,'ETAS2D8P')
        vParName = {'\mu' 'A' '\alpha' 'c' 'p' 'd' 'q' '\gamma'};
    end

    % consider earthquakes within target region R
    vCatTR   = vCat(Model.ETAS.linTR,:);        % catalogue of events in R and during [Ts, Te]
    nEQ      = length(vCatTR(:,1));
    vCumData = (1:nEQ)';
    fCumData0 = find(vCat(Model.ETAS.inR,1) < fTs,1,'last'); % find the last index for which vTdata < fTs
    if isempty(fCumData0)
        fCumData0 = 0;
    end

    vCat_m0 = vCat;
    vCat_m0(:,4) = vCat_m0(:,4) - Model.ETAS.fM0; % subtract the reference magnitude fM0
    vCumModel = zeros(nEQ,1);
    for n = 1:nEQ
        fT = vCatTR(n,1);
        vCumModel(n) = etas2d8p_cum(fT,vETASPar,vCat_m0,Model.ETAS,Bckgrnd_est);
    end

    vt = fCumData0 + vCumModel;
    vn = fCumData0 + vCumData;
    vTmod = fCumData0 + vCumModel;
    vNmod = fCumData0 + vCumModel;
    % plot cumulative rate in transformed time
    figure('Name','ETAS: Transformed Time Diagnostic','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    fGof = 0;
    plot_rate(vt,vn,vTmod,vNmod,vParName,vETASPar,vParErr,fLle,fGof,{'transformed time','event number'});
    title(Model.sTitle{1},'Interpreter','latex','FontSize',14);
    if bSaveFig
        sFig = [sFileName,'_transformed_time'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end

    vt = vCumModel;
    vn = vCumData - vCumModel;
    vTmod = vCumModel;
    vNmod = zeros(nEQ,1);
    % plot cusum
    figure('Name','ETAS: Cusum Time Diagnostic','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    fGof = 0;
    plot_rate(vt,vn,vTmod,vNmod,vParName,vETASPar,vParErr,fLle,fGof,{'cusum time','event number'});
    title(Model.sTitle{1},'Interpreter','latex','FontSize',14);
    if bSaveFig
        sFig = [sFileName,'_cusum_time'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
    
    % plot the spatial residual
    figure('Name','ETAS: Spatial Residual','Position',vPos2);
    set(gcf,'color','w'); % background color for the figure
    pcolor(Model.vX,Model.vY,mSpatRes), shading flat;
    hold on
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        pr = plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
    end
    set(gca,'Layer','top');
    cbH = colorbar('FontSize',10);
    cbLh = get(cbH,'Label');
    set(cbLh,'String','$R^\mathrm{spat}(x,y;h)$','Interpreter','latex','FontSize',12);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    %colormap(ax,"jet");
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    sTitle = {Model.sTitle{1},'Spatial residual'};
    title(sTitle,'Interpreter','latex','FontSize',14);
    hold off
    if bSaveFig
        sFig = [sFileName,'_spatial_residual'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
end

function plot_rate(vX,vY,vXmod,vModel,vParName,vETASPar,vParErr,fLle,fGof,sXYlabels)
%
    nPar = length(vETASPar);
    sposx = 0.05; % left string position
    sposy = 0.95;
    co = get(0,'DefaultAxesColorOrder');
    pl = plot(vX,vY,'.',vXmod,vModel,'-'); %
    pl(1).Color = co(2,:);
    pl(2).Color = co(1,:);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    xlabel(sXYlabels{1},'Interpreter','latex','FontSize',14);
    ylabel(sXYlabels{2},'Interpreter','latex','FontSize',14);
    % add parameter legend
    dy = 0.05;
    for n = 1:nPar
        %str = sprintf('%s = %g \xB1 %g',vParName{n},vETASPar(n),vParErr(n));
        str = ['$',vParName{n},' = ',num2str(vETASPar(n),'%.4g'),' \pm ',num2str(vParErr(n),'%.4g'),'$'];
        text(sposx,sposy-(n-1)*dy,str,'Units','normal','Interpreter','latex');
    end
    AIC = -2.0*fLle + 2*nPar; % standard definition
    str = sprintf('LL = %.4f; AIC = %.4f',fLle,AIC);
    text(sposx,sposy-nPar*dy,str,'Units','normal','Interpreter','latex');
    %fGof = etas_gof(vCat,ETAS,vETASPar,'PointProcess',sPointProc,'RateNorm',sRateNorm);
    %str = sprintf('gof = %f',fGof);
    %text(sposx,sposy-(nPar+1)*dy,str,'Units','normal');
end

