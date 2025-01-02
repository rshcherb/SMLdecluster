function tCat = inregion_events(tCat,vReg_targ,options)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 11 September 2024
%   ...
%   version: 1.1.0, 28 October 2024
%
    arguments
        tCat table
        vReg_targ double
        options.MapProjection char = []  % map projection: 'eqdcylin', 'mercator'
    end
    % Get geolocations for inside region indices
    if ~isempty(options.MapProjection)
        mstruct = defaultm('eqdcylin');
        mstruct.geoid = referenceEllipsoid('grs80','km');
        mstruct = defaultm(mstruct);
    
        [lat,lon] = projinv(mstruct,1000*tCat.Long,1000*tCat.Lat);
        %vRegTarg = readmatrix(sRegTarg);
        [vLat_targ,vLon_targ] = projinv(mstruct,1000*vReg_targ(:,1),1000*vReg_targ(:,2));
    end
    %bTarg_inds = inpolygon(tCat.Lon,tCat.Lat,lon_targ,lat_targ);

    [vLat_targ, vLon_targ] = georegion(vReg_targ);
    bTarg_inds = inpolygon(tCat.Lon,tCat.Lat,vLon_targ,vLat_targ);

    tCat = tCat(bTarg_inds,:);
end
