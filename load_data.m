%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 8 June 2021
%   ...
%   version 1.2.0, 18 October 2024
%

% extracting earthquakes from a catalogue
tCat = readtable_catalog(Model.sEqCatName,'MagMin',Model.fMmin,'MagMax',Model.fMmax,'DateStart',Model.vDateStart,'Tstart',Model.fT0,'Tend',Model.fTe,...
                         'DepthMax',Model.fDepthMax,'SeismRegion',Model.vSeismReg_all);
vCat = [tCat.Time, tCat.Lat, tCat.Lon, tCat.Mag];

% convert the geographical coordinates (lat, lon) to Cartesian (x, y) coordinates for earthquake locations and regions
vCat_km = vCat;
vCat_km(:,[3,2])= coord_projection(vCat(:,2:3),'MapProjection',Model.sMapProj); % it returns [x, y] so to make vCat = [t y x mag]
size(vCat_km)

% load geo features
vLatLim = [Model.vSeismReg_all(2),Model.vSeismReg_all(3)];
vLonLim = [Model.vSeismReg_all(4),Model.vSeismReg_all(5)];

% get fault shapes to plot
%qfaults = get_fault_shapes('BoundingBox',[vLonLim(1),vLatLim(1);vLonLim(2),vLatLim(2)],'FaultShape',Model.sFaultShape);
%qfaults = get_fault_shapes('BoundingBox',[vLonLim(1),vLatLim(1);vLonLim(2),vLatLim(2)]);

% get coastlines
sGSHHSfile = 'gshhs_i.b';
%indexfile = gshhs(sGSHHSfile,'createindex');
S = gshhs(sGSHHSfile,vLatLim,vLonLim);
levels = [S.Level];
L1 = S(levels == 1);
% get country borders
sGSHHSfile = 'wdb_borders_i.b';
%indexfile = gshhs(sGSHHSfile,'createindex');
Sb = gshhs(sGSHHSfile,vLatLim,vLonLim);
levels = [Sb.Level];
Lb = Sb(levels == 1);
cGeoFeatures.Coastline = L1;
cGeoFeatures.CountryBorders = Lb;
cGeoFeatures.TargetRegion = Model.vSeismReg_targ(2:end,:);
