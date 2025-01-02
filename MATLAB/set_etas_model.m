function Model = set_etas_model(vCat,Model,mIP,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 6 May 2021
%   ...
%   version 1.3.0, 10 June 2024
%
    %sPointProc = 'ETAS';    % 'ETASFrack', 'ETASFrackConv', 'ETAS_RS', 'ETASRemTrig', 'ETASRemTrigMOL'
    sOptMeth   = 'fmincon'; % local fmincon
    fTshift    = 0.0; % shift of the time axis to the left
    for k = 1:length(varargin)
%         if strcmp('PointProcess',varargin{k})
%             sPointProc = varargin{k+1};
%         end
        if strcmp('OptimMethod',varargin{k})
            sOptMeth = varargin{k+1};
        end
        if strcmp('Tshift',varargin{k})
            fTshift = varargin{k+1};
        end
    end
%   T0 ------- Ts ---------------- Te ----- T1
    ETAS.fT0      = Model.fT0;     % 
    ETAS.fTs      = Model.fTs;     % 
    ETAS.fTe      = Model.fTe;     %
    ETAS.fT1      = Model.fT1;     % 
    ETAS.fM0      = Model.fM0;     % the reference magnitude
    ETAS.fMc      = Model.fMc;     % the lower magnitude cutoff
    ETAS.fMmin    = Model.fMmin;   %
    ETAS.fMmax    = Model.fMmax;   %
    ETAS.mIP      = mIP;           % the initial parameters
    ETAS.sOptMeth = sOptMeth;       % optimization method
    ETAS.fTshift  = fTshift; 
    ETAS.fAlpha   = Model.fAlpha;  % confidence level

    ETAS.nJs      = find(vCat(:,1) >= Model.fTs,1,'first'); % the first event whose time is greater than or equal to fTs
    ETAS.nJe      = find(vCat(:,1) <= Model.fTe,1,'last');  % the last event whose time is less than or equal fTe
    
    ETAS.nNtot    = length(vCat(:,1));                      % number of events in the catalogue vTM
    ETAS.nNum     = ETAS.nJe - ETAS.nJs + 1;                % the number of earthquakes between fTs and fTe
    Model.fMagMean = mean(vCat(ETAS.nJs:ETAS.nJe,2));        % sample mean magnitude only for magnitudes between nJs and nJe
    
    Model.ETAS    = ETAS;
    
    Model.nJs     = ETAS.nJs;
    Model.nJe     = ETAS.nJe;
    %Model.nNum    = length(vTM(Model.nJs:Model.nJe,1)); % number of earthquakes between fTs and fTe
    Model.nNum    = ETAS.nNum; % number of earthquakes between fTs and fTe

    if ~isfield(Model,'sPointProc')
        Model.sPointProc = 'ETAS';  % the default point process
    end
    if strcmp(Model.sPointProc,'ETASFrackConv') || strcmp(Model.sPointProc,'ETAS_MultConv')
        if strcmp(Model.sConvKernel,'Pareto')
            Model.ETAS.nConvKernel = 1;
        elseif strcmp(Model.sConvKernel,'ExpPdf')
            Model.ETAS.nConvKernel = 2;
        elseif strcmp(Model.sConvKernel,'PowerLaw')
            Model.ETAS.nConvKernel = 3;
        end
        Model.ETAS.sConvKernel = Model.sConvKernel;
    end
end
