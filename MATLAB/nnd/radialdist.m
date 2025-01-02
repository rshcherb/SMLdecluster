function rdist = radialdist(lat1,lon1,lat2,lon2)
%
%   computes the radial distance in km between two points on a sphere
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.0.0, 25 May 2022
%
    fR_E = 6371.0;
    rad = pi/180.0;
%     rdist = fR_E*acos(sin(lat1*rad)*sin(lat2*rad) + cos(lat1*rad)*cos(lat2*rad)*cos((lon1-lon2)*rad));
    rdist = 2*fR_E*asin(sqrt(sin(0.5*rad*(lat1-lat2))^2 + cos(rad*lat1)*cos(rad*lat2)*sin(0.5*rad*(lon1-lon2))^2)); % Haversine formula
end