function etas2d_simulator(fBeta,vPar,fMmin,fM0,fT0,fT1,nNmax,vReg,varargin)
%
%   fBeta  - \beta = log(10)*b, b-value for the GR scaling law
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%   fMmin  - the minimum magnitude cutoff it is assumed that fMmin = fMc - fDm
%   fM0    - the reference magnitude
%   fT0    - start time of simulations
%   fT1    - the end time interval after which no more triggerings are initiated
%   nNmax  - the absolute mixumum number of events to simulate
%   vReg   - the region te generate seismicity: convex polygon, with vertices counterclockwise starting from the lower right corner
%   vEqCat - global variable for the event catalogue
%   vEqNum - global variable for the total event count
%
%        vTM            vTMout
%   fT0 ------- fTs ---------------- fTe
%                  future time window
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 April 2021
%   ...
%   version: 1.0.1, 4 January 2022
%
    global vEqCat;     % earthquake catalogue
    global nEqNum;     % earthquake counter

    sPointProc  = 'ETAS2D8P'; % the default process to simulate
    bDm         = false;
    bBackground = false;
    bAnisotropy = false;
    bSaveCat    = false;
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('Background',varargin{k})
            bBackground = true;
            vX  = varargin{k+1};
            vY  = varargin{k+2};
            mU  = varargin{k+3};
        end
        if strcmp('Anisotropy',varargin{k})
            bAnisotropy = varargin{k+1};
        end
        if strcmp('SaveCat',varargin{k})
            bSaveCat = true;
            sFileName = varargin{k+1};
        end
        if strcmp('MagBin',varargin{k})
            bDm = true;
            fDm = varargin{k+1};
        end
    end
    
    vBgCat = zeros(10000,5);
    beta1  = 1.0/fBeta;
    fEnv   = vPar(1);
    if bBackground
        fUmax = max(max(mU));   % the maximum value
        if bAnisotropy
            mAnistroProb = mU/fUmax;  % for using with etas2d_sequence_anisotropy()
        end
        fEnv = vPar(1)*fUmax;
    end
    xmin  = min(vReg(:,1)); % xmin
    Dx    = max(vReg(:,1)) - xmin; % xmax - xmin
    ymin  = min(vReg(:,2)); % ymin
    Dy    = max(vReg(:,2)) - ymin; % ymax - ymin
    vBndr = [xmin,xmin+Dx,ymin,ymin+Dy];
    
    tic
    % the thinning method: Zhuang and Touati 2012 CORSSA p. 23
    fRarea  = Dx*Dy;  % area
    fDT     = fT1 - fT0;
    fOmega  = fEnv*fRarea*fDT;
    nK      = poissrnd(fOmega);
    nBg     = 0;
    nMsRank = 1;                          % mainshock/background rank
    % generate nK background events
    for i = 1:nK
        % random position in a rectangular box
        lon = xmin + Dx*rand();           % x-coordinate
        lat = ymin + Dy*rand();           % y-coordiante
        t   = fDT*rand();
        if bBackground
            lambda = vPar(1)*interp2(vX,vY,mU,lon,lat);
            %disp(lambda/fEnv)
        else
            lambda = vPar(1);
        end
        if fEnv*rand() < lambda
            nBg = nBg + 1;
            fMag  = fMmin + exprnd(beta1);         % magnitude
            vBgCat(nBg,:) = [t,lat,lon,fMag,nMsRank];
        end
    end
    vBgCat = sortrows(vBgCat,1);
    vBgCat = vBgCat(vBgCat(:,1)>0.0,:);        % remove zero time entries

    fTmax = fT0 + fT1;
    % each background event above magnitude M0 can trigger an aftershock sequence
    for i = 1:nBg
        nEqNum = nEqNum + 1;
        vEqCat(nEqNum,:) = vBgCat(i,:);  % update the catalogue
        if vBgCat(i,4) >= fM0
            Parent.time = vBgCat(i,1);
            Parent.lat  = vBgCat(i,2);
            Parent.lon  = vBgCat(i,3);
            Parent.mag  = vBgCat(i,4);
            Parent.rank = vBgCat(i,5);
            % start the aftershock sequence
            if ~bAnisotropy
                etas2d_sequence(Parent,beta1,vPar,fMmin,fM0,fTmax,nNmax,sPointProc,vBndr); % isotropic aftershock distribution
            else
                etas2d_sequence_anisotropy(Parent,beta1,vPar,fMmin,fM0,fTmax,nNmax,sPointProc,vBndr,vX,vY,mAnistroProb); % the aftershock distribution is controlled by the background rate
            end
        end
    end
    toc

    vEqCat = vEqCat(vEqCat(:,1) > 0.0 & vEqCat(:,1) <= fT1,:);      % remove all events < 0 and > T1
    vEqCat = sortrows(vEqCat,1);               % sort the events according to increasing times
    vEqCat = vEqCat(inpolygon(vEqCat(:,3),vEqCat(:,2),vReg(:,1),vReg(:,2)),:); % keep only events within vReg
    vEqCat(:,1) = fT0 + vEqCat(:,1);
    
    if bDm   % bin magnitudes
        vEqCat(:,4) = round(vEqCat(:,4)/fDm)*fDm;
    end
    
    if bSaveCat
        fid = fopen([sFileName,'_catalog.dat'],'w');
        for n = 1:length(vEqCat(:,1))
            fprintf(fid,'%-12f  %5f  %5f  %5.3f  %4d\n',vEqCat(n,1),vEqCat(n,2),vEqCat(n,3),vEqCat(n,4),vEqCat(n,5));
        end
        fclose(fid);
    end
end



