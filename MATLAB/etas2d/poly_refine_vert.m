function [vReg_new, nVer_new] = poly_refine_vert(vReg,nNk,varargin)
%
%   Subdivides the boundary of vReg, into additional vertices 
%
%   vReg     - the original region R bounding the target events: list of vertices in counterclockwise order
%   nNk      - the initial number of radial segments, the final may differ
%
%   vReg_new - a new refined list of vetices that included the original vReg vertices
%   nVer_new - the number of new vertices subdividing the original region vReg
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 7 May 2021
%   ...
%   version: 1.0.0, 11 June 2021
%
    sMethod = 'EqualDist'; % 'EqualAngle';
    for k = 1:length(varargin)
        if strcmp('Method',varargin{k})
            sMethod = varargin{k+1};
        end
    end
    
    nVer = length(vReg(:,1)) - 1; % the last element of vReg coincides with the first
    segment = zeros(nVer,1);      % the lengths of the segments of the original polygon
    perimeter = 0;
    for k = 1:nVer
        segment(k) = sqrt( (vReg(k,1)-vReg(k+1,1))^2 + (vReg(k,2)-vReg(k+1,2))^2 ); % sqrt( (x1-x2)^2 + (y1-y2)^2 )
        perimeter = perimeter + segment(k);
    end

    if strcmp(sMethod,'EqualDist')
        % this is for events inside the target region R. The radial
        % segments are computed by connecting each event with the points on the boundary of the target region R 
        vReg_new = [];
        avd = perimeter/nNk;
        nc = 0;
        for k = 1:nVer
            nc = nc + 1;
            vReg_new(nc,1:2) = vReg(k,:);
            fvert = segment(k)/avd;
            nvert = floor(fvert);
            if (fvert - nvert) < 0.5
                nvert = nvert - 1;   % when the difference between the rounded value and initial is small don't add the last vertex
            end
            if nvert > 0
                for nv = 1:nvert
                    nc = nc + 1;
                    x = vReg(k,1) + nv*avd/segment(k)*(vReg(k+1,1)-vReg(k,1));
                    y = vReg(k,2) + nv*avd/segment(k)*(vReg(k+1,2)-vReg(k,2));
                    vReg_new(nc,:) = [x,y];
                end
            end
        end
        vReg_new(end+1,1:2) = vReg_new(1,1:2);      % add the first vertex at the end
        nVer_new = length(vReg_new(:,1)) - 1;       % do not count the last vertex
        %disp([nNk,nVer_new])
    elseif strcmp(sMethod,'EqualAngle')

    end
end
