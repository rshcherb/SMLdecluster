function [vLat, vLon] = georegion(vReg,varargin)
%
%   Computes (lat,lon) of a boundary defining the region on a map
%   vReg - the region
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 17 November 2019
%   ...
%   version 1.2.0, 1 March 2021
%
    nPts = 100;
    for k = 1:length(varargin)
        if strcmp('Points',varargin{k})
            nPts = varargin{k+1};
        end
    end
    switch vReg(1)
        case 1 % box with epicentre and size in degreees
            % (vReg(2), vReg(3) - centre; vReg(4) the size in degrees
            r = 0.5*vReg(4); % half the size
            vLat = [vReg(2)-r,vReg(2)+r,vReg(2)+r,vReg(2)-r,vReg(2)-r];
            vLon = [vReg(3)-r,vReg(3)-r,vReg(3)+r,vReg(3)+r,vReg(3)-r];
        case 2 % box (lat, lon)
            vLat = [vReg(2),vReg(3),vReg(3),vReg(2),vReg(2)];
            vLon = [vReg(4),vReg(4),vReg(5),vReg(5),vReg(4)];
        case 3 % box with the epicentre and size in km
            % (vReg(2), vReg(3)) - centre; vReg(4) - the size in km
            
        case 4 % circle
            % (vReg(2), vReg(3)) - centre; vReg(4) - the radius
            [vLat,vLon] = scircle1(vReg(2),vReg(3),vReg(4),[],[],'degrees',nPts);
        case 5 % ellipse
            % (vReg(2), vReg(3)) - centre; vReg(4) - semimajor axis; vReg(5) - semiminor axis; vReg(6) - azimuth offset
            % the semimajor axis is defined along a meridian and semiminor along a parralel
            % the azimuth offset is defined clockwise from the direction to the north
            ecc = axes2ecc(vReg(4),vReg(5));
%             [vLat,vLon] = ellipse1(vReg(2),vReg(3),[vReg(4) ecc],vReg(6));
            [vLat,vLon] = ellipse1(vReg(2),vReg(3),[vReg(4) ecc],vReg(6),[],[],'degrees',nPts);
        case 6 % polygon
            vLat = vReg(2:end,1); % the column of latitudes of polygon vertices
            vLon = vReg(2:end,2); % the column of longitudes of polygon vertices
    end
end

