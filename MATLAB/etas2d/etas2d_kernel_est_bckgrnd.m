function Bckgrnd = etas2d_kernel_est_bckgrnd(fT,vCat,Model,vPar,varargin)
%
%   Computes several measures of spatial seismicity using the variabe kernel estimate
%   (Zhuang, EPS, 2011, p.208, Eq. (10); Zhuang et al., JGR 2005, p.4, Eq. (20))
%
%   fT      - the time end of the study period
%   vCat    - event catalogue
%   Model   - model structure
%   vPar    - the ETAS model parameters
%   Bckgrnd - structure with fields:
%       vHj     - 
%       vBgProb - the probability that event i is a background/spontaneous event
%       vMu     - 
%       mMu     - variable kernel estimate of the spatial background rate
%       vBgDeclstr - 
%       mMtot   - total seismicity rate
%       mClust  - clustering coefficient
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 7 June 2021
%   ...
%   version: 1.0.0, 14 June 2021
%
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
    end
    vCat_m0 = vCat;
    vCat_m0(:,4) = vCat_m0(:,4) - Model.fM0; % subtract the reference magnitude fM0

    nJ = find(vCat_m0(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= vT(j)    
    % for each earthquake in the catalogue interpolate the baground rate vPar(1)*u(xi,yi) at each epicentre (xi, yi)
    vU = interp2(Model.vX,Model.vY,Model.mU,vCat_m0(1:nJ,3),vCat_m0(1:nJ,2));
    vBgProb_init = etas2d_bckgrnd_prob(vCat_m0,vU,vPar,fT,Model.fM0,'PointProcess',sPointProc);
    [vBgDeclstr_init, indxBg_init] = etas2d_decluster(fT,vCat,vBgProb_init);
    %
    %vHj = spat_bandwidth(vCat,nJ,Model.vBpar(1),Model.vBpar(2));
    vHj = spat_bandwidth(vCat_m0(:,[3,2]),nJ,Model.vBpar(1),Model.vBpar(2));
    [mU_est, mUtot_est, mClust_est] = kernel_est_bckgrnd_grid(Model.fT0,fT,Model.vX,Model.vY,vCat_m0,vHj,vBgProb_init);
    [vU_est, vMtot_est, vClust_est] = kernel_est_bckgrnd_events(Model.fT0,fT,vCat_m0,vHj,vBgProb_init); % estimates the background rate vUest at each earthquake location
    vBgProb_est = etas2d_bckgrnd_prob(vCat_m0,vU_est,vPar,fT,Model.fM0,'PointProcess',sPointProc);
    [vBgDeclstr_est, indxBg_est] = etas2d_decluster(fT,vCat,vBgProb_est); % one needs to provide the full catalogue with propoer magnitudes
    %
    Bckgrnd.vHj     = vHj;
    Bckgrnd.vBgProb = vBgProb_est;
    Bckgrnd.vU      = vU_est;
    Bckgrnd.mU      = mU_est;
    Bckgrnd.mMu     = vPar(1)*mU_est;
    Bckgrnd.vBgDeclstr = vBgDeclstr_est;
    Bckgrnd.indxBg  = indxBg_est;
    Bckgrnd.mUtot   = mUtot_est;
    Bckgrnd.mMutot  = vPar(1)*mUtot_est;
    Bckgrnd.mClust  = mClust_est;
    %
    Bckgrnd.vBgDeclstr_true = vBgDeclstr_init;
    Bckgrnd.indxBg_true = indxBg_init;
end
