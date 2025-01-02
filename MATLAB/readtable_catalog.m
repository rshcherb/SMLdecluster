function [tCat, indx] = readtable_catalog(sEqCatName,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 2 October 2024
%   ...
%   version 1.0.0, 2 October 2024
%
    fMmin         = 0.0;   % can be negative 
    fMmax         = 9.9;
    fDepthMin     = -10.0; % depth is negative up
    fDepthMax     = 700.0; % depth is positive down
    vDateStart    = [1900 1 1 0.0];
    vDateEnd      = [2100 1 1 0.0];
    fT0           = 0.0;
    fT1           = 365.0;
    vSeismRegion  = [2 -90 90 -360 360];
    for k = 1:length(varargin)
        if strcmp('MagMin',varargin{k})
            fMmin = varargin{k+1};
        end
        if strcmp('MagMax',varargin{k})
            fMmax = varargin{k+1};
        end
        if strcmp('DateStart',varargin{k})
            vDateStart = varargin{k+1};
        end
        if strcmp('DateEnd',varargin{k})
            vDateEnd = varargin{k+1};
        end
        if strcmp('Tstart',varargin{k})
            fT0 = varargin{k+1};
        end
        if strcmp('Tend',varargin{k})
            fT1 = varargin{k+1};
        end
        if strcmp('SeismRegion',varargin{k})
            vSeismRegion = varargin{k+1};
        end
        if strcmp('DepthMin',varargin{k})
            fDepthMin = varargin{k+1};
        end
        if strcmp('DepthMax',varargin{k})
            fDepthMax = varargin{k+1};
        end
    end
    tCat = readtable(sEqCatName);      % loading a user specified catalog file with 10 columns
    tCat.Properties.VariableNames = {'Year','Month','Day','Hour','Min','Sec','Lat','Lon','Mag','Depth'}; % if the file does not have column names
    indxm = (tCat.Mag >= fMmin) & (tCat.Mag <= fMmax) & (tCat.Depth >= fDepthMin) & (tCat.Depth <= fDepthMax);
    tCat = tCat(indxm,:);
    %
    Time = datenum(tCat.Year,tCat.Month,tCat.Day,tCat.Hour,tCat.Min,tCat.Sec) - datenum(vDateStart); % serial days relative to vDateStart
    tCat = table(Time,tCat.Lat,tCat.Lon,tCat.Mag,tCat.Depth);
    tCat.Properties.VariableNames = {'Time','Lat','Lon','Mag','Depth'}; % if the file does not have column names
    %
    indxt = (tCat.Time >= fT0)  & (tCat.Time <= fT1); % indx for earthquakes between fT0 and fT1
    tCat = tCat(indxt,:);
    %
    % compute the polygon corresponding to specific shapes
    [vLat, vLon] = georegion(vSeismRegion);
    % find points inside the polygon
    [indx_reg, on] = inpolygon(tCat.Lon,tCat.Lat,vLon,vLat);
    tCat = tCat(indx_reg,:);
    indx = [];
end

