function plot_etas2d_decluster(Backgrnd,Model,mMu,varargin)
%
%   Plots the declustered events using the backhround seismicity rate mMu
%
%   Bckgrnd    - structure for stochastic declustering
%   Model      - the model structure
%       vX     - vector for X coordinates of the renctangle
%       vY     - vector for Y coordinates of the renctangle
%       sTitle - title to add to the figure
%   mMu        - matrix for the background rate \mu*u(x,y) at the grid points (X,Y)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 16 May 2021
%   ...
%   version: 1.3.2, 18 October 2024
%
    bBackgrndRate = true;
    bDeclusteredEvents = false;
    bAftershocks = false;
    bTrueBackgrndEvents = false;
    bCoastline  = false;
    bTargetReg  = false;
    sScale      = 'Log10';  % 'Lin'
    bSaveFig    = false;
    for k = 1:length(varargin)
        if strcmp('TrueBackgrndEvents',varargin{k})
            vBgEvents_true = varargin{k+1};
            if ~isempty(vBgEvents_true)
                bTrueBackgrndEvents = true;
            end
        end
        if strcmp('DeclusteredEvents',varargin{k})
            bDeclusteredEvents = true;
            sDeclustSymbol = varargin{k+1}; % 'Dot' or 'Circle'
        end
        if strcmp('Aftershocks',varargin{k})
            bAftershocks = true;
            vAsCat = varargin{k+1};
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
        if strcmp('Scale',varargin{k})
            sScale = varargin{k+1};
        end
        if strcmp('SaveFigure',varargin{k})
            bSaveFig = true;
            sFileName = varargin{k+1};
        end
    end
    vPos = [500 300 700 600];

    if bTrueBackgrndEvents
        nBgTrue = length(vBgEvents_true(:,1));
        %disp(['Number of true background events: ',num2str(nBgTrue)])
    end
    nBgDeclstr = length(Backgrnd.vBgDeclstr(:,1));
    %disp(['Number of declustered events:     ',num2str(nBgDeclstr)])

    cbString = '$\mu\,u(x,y)$';
    if strcmp(sScale,'Log10')
        mMu = log10(mMu);
        cbString = '$\log_{10}[\mu\,u(x,y)]$';
    end
    
    figure('Name','ETAS: Estimated Background Rate and Declustered Events','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    if bBackgrndRate
        pcolor(Model.vX,Model.vY,mMu), shading flat;
        cb = colorbar('FontSize',10);
        cbLh = get(cb,'Label');
        set(cbLh,'String',cbString,'Interpreter','latex','FontSize',14);
        set(gca,'Layer','top');
        hold on
    end
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
        hold on
    end
    if bTargetReg
        plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
        hold on
    end
    if bDeclusteredEvents
        mfc = 'cyan';
        if strcmp(sDeclustSymbol,'Circle')
            vMarkerArea = (2.5 + Backgrnd.vBgDeclstr(:,4) - min(Backgrnd.vBgDeclstr(:,4))).^3;
            pb = scatter(Backgrnd.vBgDeclstr(:,3),Backgrnd.vBgDeclstr(:,2),vMarkerArea,'MarkerEdgeColor',mfc,'LineWidth',0.8);  % background events from stochastic declustering
            sFilePart = 'circle';
        elseif strcmp(sDeclustSymbol,'Dot')
            %vMarkerArea = 8;
            %pb = scatter(Backgrnd.vBgDeclstr(:,3),Backgrnd.vBgDeclstr(:,2),vMarkerArea,'filled','MarkerFaceColor',mfc,'MarkerFaceAlpha',0.7);  % background events from stochastic declustering
            pb = plot(Backgrnd.vBgDeclstr(:,3),Backgrnd.vBgDeclstr(:,2),'.','Color','cyan');  % background events from stochastic declustering
            sFilePart = 'dot';
        end
        hold on
    end
    if bAftershocks
        nAftershocks = sum(~Backgrnd.indxBg);
        pa = plot(vAsCat(:,3),vAsCat(:,2),'.','Color',"#A2142F"); % aftershocks
        hold on
    end
    if bTrueBackgrndEvents
        pbt = plot(vBgEvents_true(:,3),vBgEvents_true(:,2),'r.'); % true background events from the catalogue (ETAS simulation)
        hold on
    end
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    if bBackgrndRate && (~bDeclusteredEvents && ~bTrueBackgrndEvents && ~bAftershocks) % plot only background rate
        sFig = [sFileName,'_kernel_est_background_rate'];
        sTitle = {Model.sTitle{1},'Background rate, $\mu\, u(x,y)$'};
    end
    if bDeclusteredEvents % plot only the declustered events over background rate
        sFig = [sFileName,'_kernel_est_declustered_events_',sFilePart];
        sTitle = {[Model.sTitle{1},'; $N_\mathrm{dc} = ',num2str(nBgDeclstr),'$'],...
                  'Background rate, $\mu\, u(x,y)$, and declustered events'};
    end
    if bDeclusteredEvents && bTrueBackgrndEvents
        legend([pb, pbt],{['declustered: ',num2str(nBgDeclstr)],['true: ',num2str(nBgTrue)]});
        sFig = [sFileName,'_kernel_est_true_and_declustered_events_',sFilePart];
        sTitle = {[Model.sTitle{1},'; $N_\mathrm{dc} = ',num2str(nBgDeclstr),'$'],...
                  'Background rate, $\mu\, u(x,y)$, true and declustered events'};
    end
    if bDeclusteredEvents && bAftershocks
        legend([pa, pb],{['aftershocks: ',num2str(nAftershocks)],['declustered: ',num2str(nBgDeclstr)]});
        sFig = [sFileName,'_kernel_est_aftershocks_and_declustered_events_',sFilePart];
        sTitle = {[Model.sTitle{1},'; $N_\mathrm{dc} = ',num2str(nBgDeclstr),'$'],...
                  'Background rate, $\mu\, u(x,y)$, aftershocks and declustered events'};
    end
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    title(sTitle,'Interpreter','latex','FontSize',14);
    hold off

    if bSaveFig
        save_cf(gcf,sFig,'png','pdf','fig');
    end
end



