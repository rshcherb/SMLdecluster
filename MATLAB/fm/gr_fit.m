function [vParEst, vParErr, nNum, fMagMean] = gr_fit(vMag,fMc,fMmax,fDm,varargin)
%   [vParEst, vParErr, nNum, fMagMean] = gr_fit() - computes parameters (b-value and a) of the Gutenberg-Richter scaling relation
%
%   vMag  - earthquake magnitudes
%   fMc   - lower magnitude cutoff above which to use magnitudes
%   fMmax - upper magnitude cutoff
%   fDm   - catalog binning
%
%   Author:  Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, June 30, 2008
%   ...
%   version: 5.1.0, 2 May 2020
%
    sFitMethod = 'GuttorpHopkins'; %'Bender'; % 'Utsu'; % 
    sErrMethod = 'TintiMulargia'; % 'Aki', 'ShiBolt'
    fAlpha     = 0.05;
    sDisplay   = 'off';
    for k = 1:length(varargin)
        if strcmp('FitMethod',varargin{k})
            sFitMethod = varargin{k+1};
        end
        if strcmp('ErrorMethod',varargin{k})
            sErrMethod = varargin{k+1};
        end
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
        if strcmp('Display',varargin{k})
            sDisplay = varargin{k+1};
        end
    end
    
    indx = (vMag >= fMc) & (vMag <= fMmax); % extract all earthquakes between fMc and fMmax
    vMagDec  = vMag(indx);
    fMagMean = mean(vMagDec);      % sample mean magnitude
    nNum     = length(vMagDec);    % number of earthquakes between fMc and fMmax

    % estimate the b-value
    if strcmp(sFitMethod,'Utsu')
        fBval = b_val_utsu(fMagMean,fMc,fDm);
        sMethodTitle = '(Utsu 1965; Guo and Ogata 1997)';
    elseif strcmp(sFitMethod,'Bender')
        fBval = b_val_bender(fMagMean,fMc,fDm);
        sMethodTitle = '(Bender 1983)';
    elseif strcmp(sFitMethod,'GuttorpHopkins')
        fBval = b_val_gh(fMagMean,fMc,fDm);
        sMethodTitle = '(Guttorp and Hopkins 1986)';
    end

    % estimate the errors at a given confidence level fAlpha
    if strcmp(sErrMethod,'TintiMulargia')
        fBerr = b_err_tm(nNum,fMagMean,fMc,fDm,fAlpha);
        sErrorTitle = '(Tinti and Mulargia 1987)';
    elseif strcmp(sErrMethod,'Aki')
        fBerr = b_err_aki(nNum,fBval,fAlpha);
        sErrorTitle = '(Aki 1965)';
    elseif strcmp(sErrMethod,'ShiBolt')
        fBerr = b_err_shibolt(nNum,fBval,fAlpha);
        sErrorTitle = '(Utsu 1966; Shi and Bolt 1982)';
    end
    
    if strcmp(sDisplay,'on')
        disp(['b-value ',sMethodTitle]);
        disp(['Confidence interval ',sErrorTitle,' at ',num2str(100*(1-fAlpha)),'% confidence level']);
        disp(['    b-value = ',num2str(fBval),' +/- ',num2str(fBerr)]);
        fprintf('\n');
    end
    
    vParEst(1) = fBval;
    vParErr(1) = fBerr;
    % computing the a-value and the corresponding error
    vParEst(2) = log10(nNum) + vParEst(1)*fMc;
    vParErr(2) = vParEst(1)*fMc*sqrt((vParErr(1)/vParEst(1))^2+(fDm/fMc)^2);
end

% formula for b-value by Bender (BSSA v73, 1983, p831) for mmax = infinity
function fBval = b_val_bender(fMagMean,fMc,fDm)
%
    %z = (fMagMean - fMc)/fDm;
    %fPest = z/(z + 1.0)
    fPest = (fMagMean - fMc)/(fMagMean - fMc + fDm);
    fBeta = -log(fPest)/fDm;
    %fBval = log10(exp(1.0))*fBeta;
    fBval = 0.434294481903252*fBeta;
end

% formula for b-value by Guttorp and Hopkins (BSSA v76, 1986, p890)
function fBval = b_val_gh(fMagMean,fMc,fDm)
%
    fPest = 1 + fDm/(fMagMean - fMc);
    fBeta = log(fPest)/fDm;
    %fBval = log10(exp(1.0))*fBeta;
    fBval = 0.434294481903252*fBeta;
end

% confidence intervals from Tinti and Mulargia (BSSA v77, 1987, p2125)
function fBerr = b_err_tm(nNum,fMagMean,fMc,fDm,fAlpha)
%
    fPest = (fMagMean - fMc)/(fMagMean - fMc + fDm);
    fPCalpha = norminv(1.0-0.5*fAlpha,0,1);
    db = (1.0 - fPest)/(fDm*sqrt(nNum*fPest));
    %fBerr = log10(exp(1.0))*fPCalpha*db;
    fBerr = 0.434294481903252*fPCalpha*db;
end

% corrected formula by Utsu, 1965
% Guo and Ogata (JGR, v102, 1997, p2859) 
function fBval = b_val_utsu(fMagMean,fMc,fDm)
%
    fBeta = 1.0/(fMagMean - fMc + 0.5*fDm);
    %fBval = log10(exp(1.0))*fBeta;
    fBval = 0.434294481903252*fBeta;
end

% corrected formula by Utsu, 1965
% Guo and Ogata (JGR, v102, 1997, p2859) 
function fBerr = b_err_aki(nNum,fBval,fAlpha)
%
    fPCalpha = norminv(1.0-0.5*fAlpha,0,1);
    fBerr = fPCalpha*fBval/sqrt(nNum);
end

% Confidence intervals based on Utsu (1966) and Shi and Bolt (1982)
% Tinti and Mulargia (BSSA v77, 1987, p2127)
function fBerr = b_err_shibolt(nNum,fBval,fAlpha)
%
    fQalpha = chi2inv(1.0-0.5*fAlpha,2*nNum);
    fBerr = fBval*(fQalpha/(2*nNum) - 1.0);
end
