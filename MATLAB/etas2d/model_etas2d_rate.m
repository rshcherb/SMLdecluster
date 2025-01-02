function [vPar_est, vParErr, fLle, Backgrnd_est] = model_etas2d_rate(vCat,Model,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 31 May 2021
%   ...
%   version 2.0.0, 25 November 2023
%
    sPointProc   = 'ETAS2D8P';   % 
    sBackground  = 'Uniform';       %'Fractal';     %   'Catalog';  %
    bTrueBckgrndEvents = false;
    vBgEvents_true = [];
    sRateNorm    = 'Zhuang';     % 'Ogata', : ETAS rate normalization form
    sErrorMethod = 'Hessian';
    sSolver      = 'fmincon';    % 'gs'; % 'fminsearch'; % 
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('Background',varargin{k})
            sBackground = varargin{k+1};
        end
        if strcmp('TrueBckgrEvents',varargin{k})
            vBgEvents_true = varargin{k+1};
            if ~isempty(vBgEvents_true)
                bTrueBckgrndEvents = true;
            end
        end
        if strcmp('RateNorm',varargin{k})
            sRateNorm = varargin{k+1};
        end
        if strcmp('ErrorMethod',varargin{k})
            sErrorMethod = varargin{k+1};
        end
        if strcmp('Solver',varargin{k})
            sSolver = varargin{k+1};
        end
    end
    vCat_m0 = vCat;
    vCat_m0(:,4) = vCat_m0(:,4) - Model.ETAS.fM0; % subtract the reference magnitude fM0
    
    tic
    %[vPar_est, vParErr, fLle, exitflag, mHessian, mCovMat] = etas2d_fit(vCat,Model.ETAS,...
    %    'PointProcess',sPointProc,'RateNorm',sRateNorm,'ErrorMethod',sErrorMethod,'Solver',sSolver);
    [vPar_est, vParErr, Backgrnd_est, fLle, exitflag, mHessian] = etas2d_fit_decluster(vCat,Model.ETAS,...
        'PointProcess',sPointProc,'RateNorm',sRateNorm,'Background',sBackground,'ErrorMethod',sErrorMethod,'Solver',sSolver);
    toc
    disp(['Exit flag: ',num2str(exitflag)]);
    %save([sFileName,'_eq_cat.dat'],'vCat','-ascii');
    vParSave = [vPar_est; vParErr];
    save([Model.sFileName,'_parameters.txt'],'vParSave','-ascii');

    %Model.mMu = vPar_est(1)*Model.mU;
    fT = Model.fTe;
    if bTrueBckgrndEvents
        % stochastic declustering using the true background rate
        vBgProb_true = etas2d_bckgrnd_prob(vCat_m0,Model.ETAS.vU,vPar_est,fT,Model.fM0,'PointProcess',sPointProc);
        [vBgDeclstr_true, indxBg_true] = etas2d_decluster(fT,vCat,vBgProb_true);
%         Bckgrnd_true.vBgProb    = vBgProb_true;
%         Bckgrnd_true.vBgDeclstr = vBgDeclstr_true;
%         Bckgrnd_true.indxBg     = indxBg_true;
    end
    [vBgDeclstr_est, indxBg_est] = etas2d_decluster(fT,vCat,Backgrnd_est.vBgProb);
    [mU_est, mUtot_est, mClust_est] = kernel_est_bckgrnd_grid(Model.fT0,fT,Model.vX,Model.vY,vCat,Model.ETAS.vHj,Backgrnd_est.vBgProb);
    
    Backgrnd_est.mMu    = vPar_est(1)*mU_est;
    Backgrnd_est.vBgDeclstr = vBgDeclstr_est;
    Backgrnd_est.indxBg = indxBg_est;
    Backgrnd_est.mMutot = vPar_est(1)*mUtot_est;
    Backgrnd_est.mClust = mClust_est;

    if strcmp(Model.sMapUnit,'degree')  % if true then convert into degrees
        vCat(:,2:3) = coord_projection(vCat(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % vCat = [lat, lon] input [x, y]
        Backgrnd_est.vBgDeclstr(:,2:3) = coord_projection(Backgrnd_est.vBgDeclstr(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % Bckgrnd_est.vBgDeclstr = [lat, lon] input [x, y]
        if ~isempty(vBgEvents_true)
            vBgEvents_true(:,2:3) = coord_projection(vBgEvents_true(:,[3,2]),'MapProjection',Model.sMapProj,'Direction','inverse'); % Bckgrnd_est.vBgDeclstr = [lat, lon] input [x, y]
        end
%         latlon = coord_projection(Model.vReg_all,'MapProjection',Model.sMapProj,'Direction','inverse');
%         Model.vReg_all = [latlon(:,2), latlon(:,1)]; % for plotting 
%         latlon = coord_projection(Model.vReg_targ,'MapProjection',Model.sMapProj,'Direction','inverse');
%         Model.vReg_targ = [latlon(:,2), latlon(:,1)]; % for plotting 
        latlon = coord_projection([Model.vX, Model.vY(1)*ones(Model.nX,1)],'MapProjection',Model.sMapProj,'Direction','inverse'); % x-coordinate is changing with fixed y-coordinate
        Model.vX = latlon(:,2); % longitude
        latlon = coord_projection([Model.vX(1)*ones(Model.nY,1), Model.vY],'MapProjection',Model.sMapProj,'Direction','inverse');  % y-coordinate is changing with fixed x-coordinate
        Model.vY = latlon(:,1); % latitude
    end

    % save the data files:
    save_etas2d_estimates(vCat,Model,Backgrnd_est,'SaveResults',Model.sFileName);
end

