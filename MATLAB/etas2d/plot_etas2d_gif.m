function plot_etas2d_gif(t,vCat,Model,mMu,vPar,varargin)
%
%   Plot the ground intensity function
%
%   t      - at which time to plot
%   vCat   - event catalogue: assumes m - m_0
%   mU     - matrix for the background rate u(x,y)
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 16 May 2021
%   ...
%   version: 1.3.1, 18 October 2024
%
    sPointProc = 'ETAS2D8P'; %
    bTargetReg = false;
    bCoastline = false;
    bSaveFig   = false;
    sScale     = 'Log10';  % 'Lin'
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
        if strcmp('Scale',varargin{k})
            sScale = varargin{k+1};
        end
    end
    vPos = [500 300 700 600];

    vX = Model.vX;
    vY = Model.vY;
    
    d2    = vPar(6)^2;
    fFac  = vPar(2)*(vPar(5)-1)/vPar(4)*(vPar(7)-1)/(pi*d2); % A*(p-1)/c*(q-1)/(\pi*d^2)
    nX    = length(vX);
    nY    = length(vY);
    mRate = zeros(nY,nX);
    if strcmp(sPointProc,'ETAS2D7P')
        for iy = 1:nY
            for ix = 1:nX
                mRate(iy,ix) = gif_etas2d7p(t,vX(ix),vY(iy),mMu(iy,ix),vPar,d2,fFac,vCat);
            end
        end
    elseif strcmp(sPointProc,'ETAS2D8P')
        for iy = 1:nY
            for ix = 1:nX
                mRate(iy,ix) = gif_etas2d8p(t,vX(ix),vY(iy),mMu(iy,ix),vPar,d2,fFac,vCat);
            end
        end
    end
    
    % clip the range
    indxR = mRate > Model.fRateMax;
    mRate(indxR) = Model.fRateMax;
    
    cbString = '$\lambda(t,x,y|H)$';
    if strcmp(sScale,'Log10')
        mRate = log10(mRate);
        cbString = '$\log_{10}[\lambda(t,x,y|H)]$';
    end
    if strcmp(Model.sMapUnit,'degree')  % if true then convert into degrees
        latlon = coord_projection([vX, vY(1)*ones(nX,1)],'MapProjection',Model.sMapProj,'Direction','inverse');
        vX = latlon(:,2);
        latlon = coord_projection([vX(1)*ones(nY,1), vY],'MapProjection',Model.sMapProj,'Direction','inverse');
        vY = latlon(:,1);
    end

    % plot the figure
    figure('Name','ETAS ground intensity','Position',vPos);
    set(gcf,'color','w'); % background color for the figure
    pcolor(vX,vY,mRate), shading flat;
    hold on
    if bCoastline
        plot([cGeoFeatures.Coastline.Lon],[cGeoFeatures.Coastline.Lat],'Color','k','LineWidth',0.8);
    end
    if bTargetReg
        pr = plot(cGeoFeatures.TargetRegion(:,2),cGeoFeatures.TargetRegion(:,1),'--','Color',"#4DBEEE",'LineWidth',1.6);
    end
    set(gca,'Layer','top');
    cb = colorbar('FontSize',10);
    cbLh = get(cb,'Label');
    set(cbLh,'String',cbString,'Interpreter','latex','FontSize',14);
    ax = gca; ax.XMinorTick = 'on'; ax.YMinorTick = 'on';
    ax.TickDir = 'out';
    if strcmp(Model.sMapUnit,'degree')
        xlabel('longitude','Interpreter','latex','FontSize',14);
        ylabel('latitude','Interpreter','latex','FontSize',14);
    elseif strcmp(Model.sMapUnit,'km')
        xlabel('easting (km)','Interpreter','latex','FontSize',14);
        ylabel('northing (km)','Interpreter','latex','FontSize',14);
    end
    sTitle = Model.sTitle;
    sTitle{2} = ['Conditional rate, $\lambda(t,x,y|H)$, at $t = ',num2str(t),'$'];
    title(sTitle,'Interpreter','latex','FontSize',14);
    hold off

    if bSaveFig
        sFig = [sFileName,'_gif_',num2str(t)];
        save_cf(gcf,sFig,'png','pdf','fig');
    end
end



