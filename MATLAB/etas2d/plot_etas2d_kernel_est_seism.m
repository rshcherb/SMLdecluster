function plot_etas2d_kernel_est_seism(vCat,Model,Backgrnd,varargin)
%
%   Plots several measures of spatial seismicity (Zhuang, EPS, 2011, p.208, Eq. (10); Zhuang et al., JGR 2005, p.4, Eq. (20))
%   by using the variabe kernel estimate
%     - variable kernel estimate of the spatial background rate
%     - total seismicity rate
%     - clustering coefficient
%
%   vCat    - event catalogue
%   Bckgrnd - the structure that defines the background rate 
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 16 May 2021
%   ...
%   version: 1.2.1, 18 October 2024
%
    bCoastline = false;
    bSaveFig   = false;
    bTargetReg = false;
    sScale     = 'Log10';  % 'Lin'
    for k = 1:length(varargin)
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

    cbStringMutot = '$\lambda(t,x,y|H)$';
    mMutot = Backgrnd.mMutot;
    if strcmp(sScale,'Log10')
        mMutot = log10(mMutot);
        cbStringMutot = '$\log_{10}[\lambda(t,x,y|H)]$';
    end
    
    % plot the estimated background residual
    figure('Name','ETAS: Estimated Background Residual','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    %pcolor(Model.vX,Model.vY,Bckgrnd.mMu - Model.mMu);
    pcolor(Model.vX,Model.vY,log10(Backgrnd.mMu./Model.mMu)), shading flat;
    set(gca,'Layer','top');
    hold on
    cbh = colorbar('FontSize',10);
    %cbH.Label.String = '$\log_{10}[\mu\, u(x,y)/\mu\, u_0(x,y)]$'; %
    %cbH.Label.Interpreter = 'latex';
    cbLh = get(cbh,'Label');
    set(cbLh,'String','$\log_{10}[\mu\, u(x,y)/\mu\, u_0(x,y)]$','Interpreter','latex','FontSize',12);
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
    end
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
    sTitle = {Model.sTitle{1},'Background residual'};
    title(sTitle,'Interpreter','latex','FontSize',14);
    hold off
    if bSaveFig
        sFig = [sFileName,'_kernel_est_background_residual'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
    
    figure('Name','ETAS: Estimated Total Intensity','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    pcolor(Model.vX,Model.vY,mMutot), shading flat;
    set(gca,'Layer','top');
    hold on
    cbh = colorbar('FontSize',10);
    cbLh = get(cbh,'Label');
    set(cbLh,'String','$\hat{\Lambda}(x,y)$','Interpreter','latex','FontSize',14);
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
    end
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    title({Model.sTitle{1},'Total intensity: $\hat{\Lambda}(x,y)$'},'Interpreter','latex','FontSize',14);
    %title('Total intensity: $\hat{\Lambda}$','Interpreter','latex');
    hold off
    if bSaveFig
        sFig = [sFileName,'_kernel_est_total_intensity'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
    
    figure('Name','ETAS: Estimated Clustering Coefficient','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    pcolor(Model.vX,Model.vY,Backgrnd.mClust), shading flat;
    set(gca,'Layer','top');
    hold on
    cbh = colorbar('FontSize',10);
    %cbH.Label.String = '$\hat{\omega}(x,y)$'; %
    %cbH.Label.Interpreter = 'latex';
    cbLh = get(cbh,'Label');
    set(cbLh,'String','$\hat{\omega}(x,y)$','Interpreter','latex','FontSize',14);
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
    end
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    title({Model.sTitle{1},'Clustering coefficient: $\hat{\omega}(x,y)$'},'Interpreter','latex','FontSize',14);
    hold off
    if bSaveFig
        sFig = [sFileName,'_kernel_est_cluster_coeff'];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
end
