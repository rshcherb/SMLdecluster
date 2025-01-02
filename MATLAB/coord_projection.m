function vVect_tr = coord_projection(vVect,varargin)
%
%   Transforms between geographic coordinates (lat,lon) and Cartesian coordinates (x,y)
%
%   vVect - vector of points given as (lat,lon) or (x,y) to convert
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 9 June 2021
%   ...
%   version: 1.1.0, 18 December 2021
%
    sDirection = 'forward';  % degree -> km, 'inverse' km -> degree
    sCartUnit  = 'km';       % 'm', units for Cartesian distance
    sMapProj   = 'eqdcylin'; % 'mercator';  % 
    for k = 1:length(varargin)
        if strcmp('Direction',varargin{k})
            sDirection = varargin{k+1};
        end
        if strcmp('Unit',varargin{k})
            sCartUnit = varargin{k+1};
        end
        if strcmp('MapProjection',varargin{k})
            sMapProj = varargin{k+1};
        end
    end
    mstruct = defaultm(sMapProj);
    %mstruct.origin = [mean(vCat(:,2)) mean(vCat(:,3)) 0];
    if strcmp(sCartUnit,'km')
        mstruct.geoid = referenceEllipsoid('grs80','km');
    elseif strcmp(sCartUnit,'m')
        mstruct.geoid = referenceEllipsoid('grs80','m');
    end
    mstruct = defaultm(mstruct);
    
    vVect_tr = vVect;
    if strcmp(sDirection,'forward')
        [vVect_tr(:,1), vVect_tr(:,2)] = projfwd(mstruct,vVect(:,1),vVect(:,2)); % vVect = [lat, lon] output [x, y]
        if strcmp(sCartUnit,'km')
            vVect_tr(:,1:2) = 0.001*vVect_tr(:,1:2); % convert to km
        end
    elseif strcmp(sDirection,'inverse')
        if strcmp(sCartUnit,'km')
            vVect(:,1:2) = 1000.0*vVect(:,1:2); % convert to km
        end
        [vVect_tr(:,1),vVect_tr(:,2)] = projinv(mstruct,vVect(:,1),vVect(:,2));  % vVect_tr = [lat, lon] input [x, y]
    end
end

