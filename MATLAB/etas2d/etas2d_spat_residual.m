function mRes = etas2d_spat_residual(Model,mMuU,vCat_m0,vPar,vReg,varargin)
%
%   Computes the spatial residual for diagnostic testing (Baddeley aet al, 2015 p. 346-349)
%
%   Model.fTs, fTe - start and end time of the study interval
%   Model.vX, vY   - vectors for X abd Y coordinates of the renctangular grid
%   mMuU           - matrix for the background rate u(x,y)
%   vCat_m0        - event catalogue with magnitudes m - m_0
%   vPar           - parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 24 December 2021
%   ...
%   version: 1.1.0, 29 December 2021
%
    sPointProc   = 'ETAS2D8P'; % 
    sQuadScheme  = 'voronoi';  % 'rectangular' 
    sDummyPoints = 'boundary'; % 'meshgrid'
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('QuadScheme',varargin{k})
            sQuadScheme = varargin{k+1};
        end
        if strcmp('DummyPoints',varargin{k})
            sDummyPoints = varargin{k+1};
            nXreg = varargin{k+2};
            nYreg = varargin{k+3};
        end
    end

    xmin = min(vReg(:,1));
    xmax = max(vReg(:,1));
    ymin = min(vReg(:,2));
    ymax = max(vReg(:,2));
    vXreg = linspace(xmin,xmax,nXreg)';  % the x-coordinates of the background grid
    vYreg = linspace(ymin,ymax,nYreg)';  % the y-coordinates of the background grid
    if strcmp(sDummyPoints,'boundary')
        nP_dum = 2*(nXreg) + 2*(nYreg-2);  % number of dummy points along the boundary of vReg
        vP_dum = zeros(nP_dum,2);
        for ix = 1:nXreg
            vP_dum(ix,1) = vXreg(ix);
            vP_dum(ix,2) = ymin;               % y of the lower side of the rectangular box
            vP_dum(nXreg + ix,1) = vXreg(ix);
            vP_dum(nXreg + ix,2) = ymax;          % y of the upper side of the rectangular box
        end
        for iy = 1:nYreg-2
            vP_dum(2*nXreg + iy,1) = xmin;        % x of the left side of the rectangular box
            vP_dum(2*nXreg + iy,2) = vYreg(iy+1);    % 
            vP_dum(2*nXreg + nYreg-2 + iy,1) = xmax; % x of the right side of the rectangular box
            vP_dum(2*nXreg + nYreg-2 + iy,2) = vYreg(iy+1); 
        end
    elseif strcmp(sDummyPoints,'meshgrid')
        nP_dum = nYreg*nXreg;
        vP_dum = zeros(nP_dum,2);
        n = 1;
        for iy = 1:nYreg
            for ix = 1:nXreg
                vP_dum(n,1) = vXreg(ix);
                vP_dum(n,2) = vYreg(iy);
                n = n + 1;
            end
        end
    end
    
%     indxP = find(vCat_m0(:,1) <= fTe); % find all earthquakes defore fTe
%     vCatP = vCat_m0(indxP,:);

    indxTse = vCat_m0(:,1) >= Model.fTs & vCat_m0(:,1) <= Model.fTe; % find all earthquakes between fTs and fTe
    vCatTse = vCat_m0(indxTse,:);
    % eliminate duplicate events from vCat or randomize their locations
    [~,indxUni,~] = unique([vCatTse(:,3),vCatTse(:,2)],'rows','stable');
    % indx_orig = (1:numel(vCatTse(:,3)))';
    % indx_dupl = setdiff(indx_orig,ia);
    % vCatP(indx_dupl,2:3) = [vCatTse(indx_dupl,2) + rand(length(indx_dupl),1)*0.001, vCatTse(indx_dupl,3) + rand(length(indx_dupl),1)*0.001]; % randomize the locations of the 
    vP = vCatTse(indxUni,[3,2]);
    nP_eq = length(vP(:,1));
    vP = [vP; vP_dum]; % combine vCatP with events along the vReg
    
    if strcmp(sQuadScheme,'voronoi')
        [~,~,w] = voronoi_cell_area(vP,vReg);
    elseif strcmp(sQuadScheme,'rectangular')
        w = 0; % needs to account for multiple points in a single cell
    end
    
    nP = nP_eq + nP_dum;
    vRes = zeros(nP,1);
    lamb = zeros(nP,1);
    h    = ones(nP,1);
    for j = 1:nP_eq
        muU = interp2(Model.vX,Model.vY,mMuU,vP(j,1),vP(j,2));
        lamb(j) = gif_etas2d8p_spat(Model.fTs,Model.fTe,Model.nJs,Model.nJe,vP(j,1),vP(j,2),muU,vPar,vCat_m0);
        vRes(j) = h(j)*lamb(j)*(1 - w(j));
    end

    for j = 1:nP_dum
        muU = interp2(Model.vX,Model.vY,mMuU,vP_dum(j,1),vP_dum(j,2));
        n = nP_eq + j;
        lamb(n) = gif_etas2d8p_spat(Model.fTs,Model.fTe,Model.nJs,Model.nJe,vP_dum(j,1),vP_dum(j,2),muU,vPar,vCat_m0);
        vRes(n) = h(n)*lamb(n)*(-w(n));
    end
    
    [xi,yi] = meshgrid(Model.vX,Model.vY);
    mRes = griddata(vP(:,1),vP(:,2),vRes,xi,yi,'cubic');
    
%     vHj = spat_bandwidth(vP,nP,Model.ETAS.vBpar(1),Model.ETAS.vBpar(2));
%     mRes = zeros(Model.nY,Model.nX);
%     for iy = 1:Model.nY
%         for ix = 1:Model.nX
%             for j = nP
%                 r = sqrt((Model.vX(ix) - vP(j,1))^2 + (Model.vY(iy) - vP(j,2))^2);
%                 fZj = gauss_kernel(r,vHj(j));
%                 mRes(iy,ix) = mRes(iy,ix) + fZj*vRes(j);
%             end
%         end
%     end
end

