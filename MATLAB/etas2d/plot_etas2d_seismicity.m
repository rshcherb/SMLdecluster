function plot_etas2d_seismicity(vCat,Model,varargin)
%
%   Plot seismicity, faults, strong earthquakes, target region, and initial background rate Model.mMu
%   vCat - event catalogue
%   Model - the model structure
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 28 April 2021
%   ...
%   version: 1.3.1, 18 October 2024
%
    bSaveFig    = false;
    bTrueBckgrndEvents = false;
    bMs         = false;
    fMsMin      = 5.0;
    bFaults     = false;
    bCoastline  = false;
    bTargetReg  = false;
    sScale      = 'Log10';  % 'Lin'
    for k = 1:length(varargin)
        if strcmp('SaveFigure',varargin{k})
            bSaveFig = true;
        end
        if strcmp('TrueBckgrEvents',varargin{k})
            vBgEvents_true = varargin{k+1};
            if ~isempty(vBgEvents_true)
                bTrueBckgrndEvents = true;
            end
        end
        if strcmp('MsMagMin',varargin{k})
            bMs = true;
            fMsMin = varargin{k+1};
        end
        if strcmp('GeoFeatures',varargin{k})
            cGeoFeatures = varargin{k+1};
            if isfield(cGeoFeatures,'TargetRegion')
                bTargetReg = true;
            end
            if isfield(cGeoFeatures,'Coastline')
                bCoastline = true;
            end
            if isfield(cGeoFeatures,'Faults')
                %bFaults = true;
            end
        end
        if strcmp('Scale',varargin{k})
            sScale = varargin{k+1};
        end
    end
    vPos = [500 300 700 600];

    cbString = '$\mu\, u_0(x,y)$';
    mMu = Model.mMu;
    if strcmp(sScale,'Log10')
        mMu = log10(Model.mMu);
        cbString = '$\log_{10}[\mu\, u_0(x,y)]$';
    end
    
    figure('Name','Seismicity Map','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    pcolor(Model.vX,Model.vY,mMu), shading flat;
    set(gca,'Layer','top');
    cb = colorbar('FontSize',10);
    cbLh = get(cb,'Label');
    set(cbLh,'String',cbString,'Interpreter','latex','FontSize',14);
    hold on
    pae = plot(vCat(:,3),vCat(:,2),'.','Color',"cyan"); % plot all earthquakes
    vpl = pae;
    clbl = {'earthquakes'};
    n = 2;
    if bTrueBckgrndEvents
        pbt = plot(vBgEvents_true(:,3),vBgEvents_true(:,2),'r.'); % true background events
        vpl = [vpl, pbt];
        clbl{n} = 'true bg';
        n = n + 1;
    end
    if bFaults
        pf = plot([vFlist(:,1)';vFlist(:,3)'],[vFlist(:,2)';vFlist(:,4)'],'Color','#A2142F','LineWidth',1.5);
        vpl = [vpl, pf(1)];
        clbl{n} = 'faults';
        n = n + 1;
    end
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        pr = plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
%         vpl = [vpl, pr];
%         clbl{n} = 'target';
%         n = n + 1;
    end
    vMarkerArea = (2.5 + vCat(:,4) - min(vCat(:,4))).^3;
    %scatter(vCat(~indxBg,3),vCat(~indxBg,2),vMarkerArea(~indxBg),'MarkerEdgeColor','black','LineWidth',0.5);  % aftershocks
    if bMs
        indxMs = vCat(:,4) >= fMsMin;            % strong events
        pms = scatter(vCat(indxMs,3),vCat(indxMs,2),vMarkerArea(indxMs),'MarkerEdgeColor',"#7E2F8E",'LineWidth',1);  % large earthquakes
        vpl = [vpl, pms];
        clbl{n} = ['mainshocks $\ge ',num2str(fMsMin,2),'$'];
        %n = n + 1;
    end
    legend(vpl,clbl,'Interpreter','latex','FontSize',10);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    title({Model.sTitle{1},'Seismicity map and initial background rate, $\mu\, u_0(x,y)$'},'Interpreter','latex','FontSize',14);
    hold off

    if bSaveFig
        sFig = [Model.sFileName,'_seismicity_map'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
end



