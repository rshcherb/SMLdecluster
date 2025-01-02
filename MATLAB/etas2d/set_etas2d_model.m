function Model = set_etas2d_model(vCat,Model,mIP,options)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 May 2021
%   ...
%   version: 1.1.0, 18 October 2024
%
    arguments
        vCat double
        Model struct
        mIP struct
        options.Solver char = 'fmincon'
        options.NumSegments int32 = 1000
    end

    %   T0 ------- Ts ---------------- Te ----- T1
    ETAS.fT0      = Model.fT0;     % 
    ETAS.fTs      = Model.fTs;     % 
    ETAS.fTe      = Model.fTe;     %
    ETAS.fT1      = Model.fT1;     % 
    ETAS.fM0      = Model.fM0;     % the reference magnitude
    ETAS.fMc      = Model.fMc;     % the lower magnitude cutoff
    ETAS.vX       = Model.vX;
    ETAS.vY       = Model.vY;
    ETAS.mU       = Model.mU;      % the initial background rate
%     ETAS.fUint    = Model.fUint;
%     ETAS.fUmean   = Model.fUmean;
%     ETAS.fUmax    = Model.fUmax;
    ETAS.vBpar    = Model.vBpar;   % the parameters for the spatial bandwith model
    ETAS.mIP      = mIP;           % the initial parameters
    ETAS.sSolver  = options.Solver; % optimization method
    ETAS.fAlpha   = Model.fAlpha;  % confidence level
    ETAS.eps1     = Model.eps1;
    ETAS.eps2     = Model.eps2;
    ETAS.eps3     = Model.eps3;

    inR           = inpolygon(vCat(:,3),vCat(:,2),Model.vReg_targ(:,1),Model.vReg_targ(:,2)); % logical indeces 
    ETAS.inR      = inR;           % logical indices of events in the target region R during [T0, Te]
    ETAS.outR     = ~inR;          % logical indeces of events outside the target region
    ETAS.linR     = find(inR);     % vector containing linear indices of each event in the target region R and during [T0,Te]
    inTR          = (vCat(:,1) >= Model.fTs & inpolygon(vCat(:,3),vCat(:,2),Model.vReg_targ(:,1),Model.vReg_targ(:,2)));
    ETAS.inTR     = inTR;          % logical indices of each event during [Ts, Te] and the target region R
    ETAS.outTR    = ~inTR;         % logical indeces of events outside the target region
    ETAS.linTR    = find(inTR);    % vector containing linear indices of each event in the target region R and during [T_s,T_e]
    ETAS.vReg_targ = Model.vReg_targ;
    
    ETAS.nJs      = find(vCat(:,1) >= Model.fTs,1,'first'); % the first event whose time is greater than or equal to fTs
    ETAS.nJe      = find(vCat(:,1) <= Model.fTe,1,'last');  % the last event whose time is less than or equal fTe
    
    ETAS.nNtot    = length(vCat(:,1));                      % number of events in the catalogue vCat
    ETAS.nNum     = ETAS.nJe - ETAS.nJs + 1;                % number of events between fTs and fTe
    ETAS.nNtarg   = length(ETAS.linTR);                     % number of events in the target region and in [Ts, Te]

    % for each event in the region create radial segments with given radius and angle
    [vRk, vRk2, vDk, vReg_ref, nVer_ref] = segment_rad_angle(Model.vReg_targ,inR,[vCat(:,3),vCat(:,2)],options.NumSegments);
    ETAS.nNk      = nVer_ref;
    ETAS.vRk      = vRk;                                    % the radius for each radial segment
    ETAS.vRk2     = vRk2;                                   % the radius for each radial segment
    ETAS.vDk      = vDk;                                    % the angle between two radial sides of a segment
    %disp([ETAS.nNk,length(ETAS.vRk),length(ETAS.vRk2),length(ETAS.vDk)])
    % for each earthquake in the catalogue interpolate the baground rate u(xi,yi) at each epicentre (xi, yi)
    ETAS.vU = zeros(ETAS.nJe,1);
    for i = 1:ETAS.nJe
        ETAS.vU(i) = interp2(ETAS.vX,ETAS.vY,ETAS.mU,vCat(i,3),vCat(i,2));
    end
    
    ETAS.vHj = spat_bandwidth(vCat(:,[3,2]),ETAS.nJe,ETAS.vBpar(1),ETAS.vBpar(2));

    % compute the volume integral of the rate function u(x,y) over the target region R
    %[ETAS.fUquadR, DT, vXp, vYp] = quad_polygon(ETAS.vX,ETAS.vY,ETAS.mU,Model.vReg_targ,'Method','Triangulation');
    vBgProb_est = ones(ETAS.nJe,1);
    ETAS.fUquadR = quad_background(ETAS.inR,ETAS.nJe,vBgProb_est,ETAS.vHj,ETAS.vRk,ETAS.vRk2,ETAS.vDk,ETAS.nNk)/(ETAS.fTe - ETAS.fT0);

    Model.ETAS    = ETAS;
    
    Model.nJs     = ETAS.nJs;
    Model.nJe     = ETAS.nJe;
    Model.nNum    = ETAS.nNum; % number of earthquakes between fTs and fTe
end


