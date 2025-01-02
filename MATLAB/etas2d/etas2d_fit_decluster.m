function [vPar_est, vParErr, Bckgrnd_est, fLle, exitflag, mHessian] = etas2d_fit_decluster(vCat,ETAS,varargin)
%
%   Fit the ETAS model to the catalogue of events using the interative stochastic declustering
%   (Zhuang, EPS, 2011, p.208, Eq. (10); Zhuang et al., JGR 2005, p.4, Eq. (20))
%
%   vCat - the catalogue of events
%   ETAS - input model structure: fT0, fTs, fTe, fT1, fM0, vU, ...
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 29 April 2021
%   ...
%   version: 1.0.2, 18 November 2023
%
    sPointProc   = 'ETAS2D8P'; % 
    sBackground  = 'None';     % 'Fractal';     %   'Catalog';  %
    sRateNorm    = 'Zhuang';   % 'Ogata'; % : ETAS rate normalization form
    sParamError  = 'on';       % whether to compute parameter errors or not: off|on
    sErrorMethod = 'Hessian';  % 'LogLik'
    sSolver      = 'fmincon';  % 'gs'; % 'fminsearch'; % 
    sDisplay     = 'off';      % 'final';    % 'iter';     % % fmincon level of display
    nIterMax     = 50;
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
        if strcmp('RateNorm',varargin{k})
            sRateNorm = varargin{k+1};
        end
        if strcmp('Background',varargin{k})
            sBackground = varargin{k+1};
        end
        if strcmp('IterMax',varargin{k})
            nIterMax = varargin{k+1};
        end
        if strcmp('ParamError',varargin{k})
            sParamError = varargin{k+1};
        end
        if strcmp('ErrorMethod',varargin{k})
            sErrorMethod = varargin{k+1};
        end
        if strcmp('Solver',varargin{k})
            sSolver = varargin{k+1};
        end
        if strcmp('Display',varargin{k})
            sDisplay = varargin{k+1};
        end
    end
    
    vCat_m0 = vCat;
    vCat_m0(:,4) = vCat_m0(:,4) - ETAS.fM0; % subtract the reference magnitude fM0
    
    vUest_prev   = ETAS.vU;
    %vMuest_prev  = ETAS.mIP.vInitParams(1)*vU; % the initial background rate
    vParEst_prev = ETAS.mIP.vInitParams;
    fLle_prev    = 1.0;
    k = 0; e1 = 1.0; e2 = 1.0; e3 = 1.0;
    while (e1 >= ETAS.eps1 || e2 >= ETAS.eps2 || e3 >= ETAS.eps3) && (k <= nIterMax)
        k = k + 1;
        if strcmp(sPointProc,'ETAS2D7P')
            if strcmp(sRateNorm,'Zhuang')
                llfun = @(vPar) -etas2d7p_llfun(vPar,vCat_m0,ETAS);
            elseif strcmp(sRateNorm,'Ogata')
            end
        elseif strcmp(sPointProc,'ETAS2D8P')
            if strcmp(sRateNorm,'Zhuang')
                llfun = @(vPar) -etas2d8p_llfun(vPar,vCat_m0,ETAS);
            elseif strcmp(sRateNorm,'Ogata')
            end
        end
        %
        [vPar_est,fLle,exitflag,mHessian] = fmin_iteration(llfun,ETAS.mIP,sSolver,sDisplay);
        %
        vBgProb_est = etas2d_bckgrnd_prob(vCat_m0,vUest_prev,vPar_est,ETAS.fTe,'PointProcess',sPointProc);
        [vU_est, vUtot, vClust] = kernel_est_bckgrnd_events(ETAS.fT0,ETAS.fTe,vCat_m0,ETAS.vHj,vBgProb_est); % estimates the background rate vUest at each earthquake location
        ETAS.vU  = vU_est;
        ETAS.mIP.vInitParams = vPar_est;
        % compute the volume integral of the rate function u(x,y) over the target region R
        if strcmp(sBackground,'Catalog')
            %[mU_est, mUtot_est, mClust_est] = kernel_est_bckgrnd_grid(ETAS.fT0,ETAS.fTe,ETAS.vX,ETAS.vY,vCat,ETAS.vHj,vBgProb_est);
            %[ETAS.fUquadR, DT, vXp, vYp] = quad_polygon(ETAS.vX,ETAS.vY,mU_est,ETAS.vReg_targ,'Method','Triangulation');
            ETAS.fUquadR = quad_background(ETAS.inR,ETAS.nJe,vBgProb_est,ETAS.vHj,ETAS.vRk,ETAS.vRk2,ETAS.vDk,ETAS.nNk)/(ETAS.fTe - ETAS.fT0);
        end
        
        e1 = max(abs(vPar_est./vParEst_prev - 1.0)); % the relative change in parameters
        %e2 = max(abs(vMu_est./vMuest_prev - 1.0));   % the relative change in the background rate
        e2 = max(abs(vU_est./vUest_prev - 1.0));     % the relative change in the background rate
        e3 = abs(fLle/fLle_prev - 1.0);              % the relative change in the log-likelihood value
        vParEst_prev = vPar_est;
        vUest_prev   = vU_est;
        fLle_prev    = fLle;
        disp(['Iteration: ',num2str(k)])
        disp(['Estimated parameters: ',num2str(vPar_est)])
        disp(['Relative change: parameters = ',num2str(e1),'; background = ',num2str(e2),'; loglik = ',num2str(e3)])
    end

    Bckgrnd_est.vBgProb = vBgProb_est;
    Bckgrnd_est.vU      = vU_est;
    Bckgrnd_est.vMu     = vPar_est(1)*vU_est;
    Bckgrnd_est.vMutot  = vPar_est(1)*vUtot;
    Bckgrnd_est.vClust  = vClust;
    
    fLle = -fLle; % change sign to reflect that the maximum of the log-likelihood function is needed

    vParErr = zeros(1,length(ETAS.mIP.vInitParams));
    % computes parameter errors for the ETAS model using one of the methods: 'Hesian', 'OVP', 'LogLik'
    mCovMat = 0;
    if strcmp(sParamError,'on')
        [vParErr, mCovMat] = etas2d_par_err(vPar_est,vCat,ETAS,mHessian,ETAS.fAlpha,...
            'ErrorMethod',sErrorMethod,'RateNorm',sRateNorm);
    end
end

function [vParEst,fLle,exitflag,mHessian] = fmin_iteration(llfun,mIP,sSolver,sDisplay)
%
    mHessian = zeros(length(mIP.vInitParams));
    if strcmp(sSolver,'fmincon')
        % problem setup for fmincon uses optimoptions() from Global search toolbox
        %options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-7,'Display',sDisplay);
        %options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-7,'Display',sDisplay,'Diagnostics','on','DiffMinChange',0.01);
        options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',1000,'MaxFunctionEvaluations',1500,'StepTolerance',1e-8,'Display',sDisplay);
        %options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',1500,'Display',sDisplay);
        problem = createOptimProblem('fmincon','objective',llfun,'x0',mIP.vInitParams,...
                            'lb',mIP.vLowB,'ub',mIP.vUppB,'Aeq',mIP.Aeq,'beq',mIP.beq,'options',options);

        [vParEst,fLle,exitflag,output,lambda,grad,mHessian] = fmincon(problem);
    elseif strcmp(sSolver,'fminsearch')
        options = optimset('Display',sDisplay,'MaxFunEvals',2000,'TolX',1e-7,'TolFun',1e-7);
        [vParEst,fLle,exitflag,output] = fminsearch(llfun,mIP.vInitParams,options);
        sParamError  = 'off';
    elseif strcmp(sSolver,'gs')
        % problem setup for global search
        options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-7,'Display',sDisplay);
        problem = createOptimProblem('fmincon','objective',llfun,'x0',mIP.vInitParams,...
                            'lb',mIP.vLowB,'ub',mIP.vUppB,'Aeq',mIP.Aeq,'beq',mIP.beq,'options',options);

        gs = GlobalSearch;
        [vParEst,fLle] = run(gs,problem);
        sParamError  = 'off';
    end

end

function [vParErr, mCovMat] = etas2d_par_err(vPar,vCat,ETAS,mHessian,fAlpha,varargin)
%
    sErrorMethod = 'Hessian';
    sRateNorm    = 'Zhuang'; % 'Ogata'; % ETAS rate normalization form
    for k = 1:length(varargin)
        if strcmp('ErrorMethod',varargin{k})
            sErrorMethod = varargin{k+1};
        end
        if strcmp('RateNorm',varargin{k})
            sRateNorm = varargin{k+1};
        end
    end
    if strcmp(sErrorMethod,'Hessian') % using the numerical Hessian from the optimization
        mFI = mHessian;
    elseif strcmp(sErrorMethod,'LogLik')
        if strcmp(sRateNorm,'Zhuang')
            %mFI = fim_etas_loglik(vPar,vTM,fTs,fTe);
        elseif strcmp(sRateNorm,'Ogata')
        %elseif strcmp(sRateNorm,'Zhuang')
        end
    end
    mCovMat = inv(mFI); % inverse of the Fisher information matrix
    fPCalpha = norminv(1-0.5*fAlpha,0,1); % (1-\alpha/2)*100% percentile of the standard normal distribution, N(0,1)

    vParErr = zeros(1,length(vPar));
    for n = 1:length(vPar)
        vParErr(n) = fPCalpha*sqrt(mCovMat(n,n));
    end
end

