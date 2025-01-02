%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 6 May 2021
%   ...
%   version 1.2.0, 11 October 2024
%
Model.sSeismName  = 'Italy';
sDir_data = [sDirPart,'catalog/'];
sDir_out  = [sDirPart,'etas_fit/'];
Model.sEqCatName  = sprintf('%sraw_%s_catalog.dat',sDir_data,Model.sRegion);
Model.vSeismReg_all = [2 35.0 48.5 5.0 21.0];
Model.vSeismReg_targ = [6 0; 36.0, 19.0; 48.0, 19.0; 48.0, 6.1; 39.5, 6.1; 36.0, 13.0; 36.0, 19.0];   % Italy region counterclockwise starting from the lower right corner
%   T0 ------- fTs ---------------- fTe ------- fT1
Model.vDateStart  = [2000,1,1,0.0,0,0]; % vector for the calendar date (mainshock) corresponding to fT0
Model.vDateS      = [2004,1,1,0.0,0,0]; % vector for the calendar date (mainshock) corresponding to fTs
Model.vDateEnd    = [2021,1,1,0.0,0,0]; % vector for the calendar date (mainshock) corresponding to fTe
Model.fT0         = 0.0;
Model.fTs         = datenum(Model.vDateS)   - datenum(Model.vDateStart); % end of the target time interval;
Model.fTe         = datenum(Model.vDateEnd) - datenum(Model.vDateStart); % end of the target time interval;
Model.fDT         = 1.0;                % forecasting time interval in days  
Model.sMapUnit    = 'degree';           % coordinates of plots in 'km' or 'degree'
Model.sMapProj    = 'eqdcylin';         % 'eqdcylin'; % 'mercator'; % map projection to use: 'mercator', 'eqdcylin'
Model.sRateNorm   = 'Zhuang'; % 'Ogata'; %  
Model.sErrorMethod= 'Hessian'; % 
sSolver           = 'fmincon'; % 'gs'; % 'fminsearch'; % 
PPip.Aeq = []; PPip.beq = [];     % no constraints on parameters
if ~exist('sPointProc','var')
    sPointProc    = 'ETAS2D8P';  % point process to use
end
Model.sPointProc  = sPointProc;
if strcmp(sPointProc,'ETAS2D7P')
    % set the initial parameters for the ETAS fitting
    %                   mu      A     alpha  c      p      d     q
    PPip.vInitParams = [1.0     0.25  1.8    0.1    1.5    0.5   1.9];
    PPip.vLowB    =    [0.0     0.0   0.0    0.0    1+1e-6 1e-6  1+1e-6];
    PPip.vUppB    =    [100.0   10.0  10.0   1.0    10.0   10.0  10.0];
elseif strcmp(sPointProc,'ETAS2D8P')
    % to introduce constraints on the parameters
    %  PPip.Aeq = zeros(8); PPip.Aeq(1,1) = 1; %
    %  PPip.beq = zeros(8,1); PPip.beq(1) = 0.0; % mu
    % set the initial parameters for the ETAS fitting
    %                   mu      A     alpha  c      p      d      q      gamma
    PPip.vInitParams = [0.5     0.25  2.0    0.05   1.2    0.5    1.5    0.1];
    PPip.vLowB    =    [0.0     0.0   0.0    0.0    1+1e-6 1e-6   1+1e-6 0.0];
    PPip.vUppB    =    [100.0   10.0  10.0   1.0    10.0   1000.0 10.0   10.0];
end
Model.sTitle0     = [replace(Model.sSeismName,"_"," "),': '];
Model.fT1         = Model.fTe;
Model.fDm         = 0.1;
Model.fAlpha      = 0.05;
Model.fMmin       = Model.fMc - 0.5*Model.fDm;  % lower cutoff for magnitudes to extract or generate
Model.fMmax       = 9.0;                        % upper cutoff for magnitudes to use
Model.fM0         = Model.fMc;                  % reference magnitude for the ETAS model
Model.fDepthMin   = -10.0;
Model.fDepthMax   = 30.0;
Model.eps1        = 0.01;
Model.eps2        = 0.1;
Model.eps3        = 0.01;
Model.nBins       = 100;          % number of bins for the magnitude x axis for computing the distributions
Model.vBinLimits  = [0,1000];     % 
Model.nMax        = 20000;        % the maximum number of events to generate

lat = Model.vSeismReg_all(2:3);
lon = Model.vSeismReg_all(4:5);
vReg_all = [lat(2) lon(1); lat(2) lon(2); lat(1) lon(2); lat(1) lon(1); lat(2) lon(1)]; % rectangular region is converted into a polygon representation
Model.vReg_all  = coord_projection(vReg_all,'MapProjection',Model.sMapProj);
Model.vReg_targ = coord_projection(Model.vSeismReg_targ(2:end,:),'MapProjection',Model.sMapProj);

bAnisotropy  = false;      % whether to use an anisotropic distribution of aftershocks 
fNormFactor  = 0.9e-6;     % the maximum value for the background rate u(x,y); the total background rate is vPar(1)*u(x,y)
Model.fH     = 0.7;        % the Hurst exponent for the fractal rate
vSmoothPar   = 1.0;        % the parameter for the Gaussian filter to smooth the background rate
if exist('sBackground','var') && strcmp(sBackground,'Catalog')
    vSmoothPar = [];   % do not smooth the background when estimating from a catlogue
end
Model.vBpar    = [5.0, 10];  % the bandwidth parameters for kernel estimates: [h_min, np]
Model.nFaults  = 15;         % the number of faults to simulate
Model.vFpar    = [1.5, 70, 0.7, 0.2]; % fault parameters: [q, d, P, I0]
Model.fRateMax = 1e-4;       % the maximum value to clip the ground intensity function
Model.nX       = 256;                        % the number of bins along x-side of the region 
Model.nY       = Model.nX;                   % the number of bins along y-axis of teh region
