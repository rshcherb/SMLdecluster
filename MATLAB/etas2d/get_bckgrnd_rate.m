function Bckgrnd = get_bckgrnd_rate(varargin)
%
%   Reads from a file or generates/estimates background seismicity rate
%
%   get_bckgrnd_rate('LoadFile',sFileName) - loading from a file
%   get_bckgrnd_rate('Uniform',Model,fNormFactor) - genrates new synthetic uniform rate
%   get_bckgrnd_rate('Fractal',Model,fNormFactor) - genrates new synthetic rate using fractional Brownian noise
%   get_bckgrnd_rate('Catalog',vCat,fT0,fT1,vReg) - estimates the background rate from a catalogue
%
%   For a new synthetic rate:
%   nXsize, nYsize - the size of the grid, for example (256,256)
%   fH             - the Hurst exponent: (0,1)
%
%   Output: Bckgrnd structure with fields:
%      vX     - vector for X coordinates of the renctangle
%      vY     - vector for Y coordinates of the renctangle
%      mU     - matrix for the background rate u(x,y)
%      vFlist - the list of faults: [vx1 vx2 vy1 vy2]
%      fQuad  - the integral of the background rate u(x,y) over the rectangular area
%      fUmean - fUint/area
%      fUmax  - max(max(mU))
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 11 May 2021
%   ...
%   version: 1.2.0, 3 January 2022
%
    bFaults   = false;
    bSmooth   = false;
    bSaveFile = false;
    sSeism    = 'All';
    for k = 1:length(varargin)
        if strcmp('LoadFile',varargin{k}) % load the background rate from a file
            sBackground = varargin{k};
            sFileName   = varargin{k+1};
        end
        if strcmp('Fractal',varargin{k})  % generate the background rate as a fractal (fractional Brownian noise)
            sBackground = varargin{k};
            Model       = varargin{k+1};
            fNormFactor = varargin{k+2};
        end
        if strcmp('Uniform',varargin{k})  % generate the uniform background rate
            sBackground = varargin{k};
            Model       = varargin{k+1};
            fNormFactor = varargin{k+2};
        end
        if strcmp('Faults',varargin{k})   % combine the fractal rate with randomly oriented faults
            sFaultsFile = varargin{k+1};
            if ~isempty(sFaultsFile)
                bFaults = true;
            end
        end
        if strcmp('Catalog',varargin{k})  % estimate the background rate from a catalogue
            sBackground = varargin{k};
            vCat        = varargin{k+1};  % the length of the catalog is controlled outside this function
            Model       = varargin{k+2};
            sSeism      = varargin{k+3};  % 'All' - all seismicity; 'NND' - declsutered; 
        end
        if strcmp('Smooth',varargin{k})   % smooh the baground rate using a Gaussian filter
            bSmooth = true;
            sMethod = varargin{k+1};
            vSmoothPar = varargin{k+2};
            if isempty(vSmoothPar)
                bSmooth = false;
            end
        end
        if strcmp('SaveFile',varargin{k})
            bSaveFile = true;
            sFileOut = varargin{k+1};
        end
    end

    vFlist = [];
    vU = [];
    if ~strcmp(sBackground,'LoadFile')
        vX = Model.vX;
        vY = Model.vY;
    end
    
    if strcmp(sBackground,'LoadFile')
        [vX, vY, mU] = read_background_rectangle(sFileName);
        if bFaults
            vFlist = read_faults(sFaultsFile);
        end
        
    elseif strcmp(sBackground,'Uniform')
        mU = fNormFactor*ones(Model.nY,Model.nX);
    elseif strcmp(sBackground,'Fractal')
        if bFaults
            [mU, vFlist] = faults_intensity(Model.vX,Model.vY,Model.nFaults,Model.vFpar,Model.fH,'NormFactor',fNormFactor);
        else
            mU = fractal_intensity(Model.nX,Model.nY,Model.fH,'NormFactor',fNormFactor);  % generate the intensity function as a fractal
        end
        if bSmooth
            if strcmp(sMethod,'GaussFilter')
                mU = imgaussfilt(mU,vSmoothPar(1));
            end
        end
        if bSaveFile
            save_background_file(Model.nX,Model.nY,Model.vX,Model.vY,mU,sFileOut);
            if bFaults
                save_faults_file(vFlist,sFaultsFile);
            end
        end
    elseif strcmp(sBackground,'Catalog')
        if strcmp(sSeism,'All')
            nJ = length(vCat(:,1));
            fT = vCat(end,1);
        elseif strcmp(sSeism,'NND')
            
        end
        
        vHj = spat_bandwidth(vCat(:,[3,2]),nJ,Model.vBpar(1),Model.vBpar(2));
        vBgProb = ones(nJ,1);     % assumes all events are background
        [mU, mMtot, mClust] = kernel_est_bckgrnd_grid(Model.fT0,fT,Model.vX,Model.vY,vCat,vHj,vBgProb);
        [vU, vMtot, vClust] = kernel_est_bckgrnd_events(Model.fT0,fT,vCat,vHj,vBgProb); % estimates the background rate vUest at each earthquake location
        if bSmooth
            if strcmp(sMethod,'GaussFilter')
                mU = imgaussfilt(mU,vSmoothPar(1));
            end
        end
    end
    fRarea = (max(vX)-min(vX))*(max(vY)-min(vY));  % area
    
    Bckgrnd.vX     = vX;
    Bckgrnd.vY     = vY;
    Bckgrnd.mU     = mU;
    Bckgrnd.vU     = vU;
    Bckgrnd.vFlist = vFlist;
    Bckgrnd.fQuad  = trapz(vY,trapz(vX,mU,2));     % integral of the intensity over the whole area
    Bckgrnd.fUmean = Bckgrnd.fQuad/fRarea;         % the mean
    Bckgrnd.fUmax  = max(max(mU));                 % the maximum value
end

% read the background rate from a file
function [vX, vY, mU] = read_background_rectangle(sFileName)
%
    fid = fopen(sFileName,'r');
    R = fscanf(fid,'%d %d %f %f %e',[5 inf]);
    R = R';                % transpose the array to match the original file
    % the input file can be in two formats: 
    %    - columns 3 and 4 are (x, y) Cortesian coordinates
    %    - columns 3 and 4 are (lon, lat) 
    nXsize = max(R(:,1));
    nYsize = max(R(:,2));
    mU = zeros(nYsize,nXsize);
    vX = zeros(nXsize,1);
    vY = zeros(nYsize,1);
    for n = 1:length(R)
        mU(R(n,2),R(n,1)) = R(n,5);
    end
    for n = 1:nXsize
        vX(n) = R(n,3);
    end
    for n = 1:nYsize
        nn = 1 + (n-1)*nXsize;
        vY(n) = R(nn,4);
    end
    %disp([vX, vY])
    fclose(fid);
end

% save the background rate into a file
function save_background_file(nX,nY,vX,vY,mU,sFileOut)
%
    fid = fopen(sFileOut,'w');
    for iy = 1:nY
        for ix = 1:nX
            fprintf(fid,'%d %d %f %f %e\n',ix,iy,vX(ix),vY(iy),mU(iy,ix));
        end
    end
    fclose(fid);
end

% save the background rate into a file
function vFaults = read_faults(sFaultsFile)
%
    fid = fopen(sFaultsFile,'r');
    F = fscanf(fid,'%f %f %f %f',[4 Inf]);
    vFaults = F';
    %disp(vFaults)
    fclose(fid);
end

% save the background rate into a file
function save_faults_file(vFaults,sFaultsFile)
%
    fid = fopen(sFaultsFile,'w');
    nF = length(vFaults(:,1));
    for n = 1:nF
        fprintf(fid,'%f %f %f %f\n',vFaults(n,1),vFaults(n,2),vFaults(n,3),vFaults(n,4));
    end
    fclose(fid);
end
