function mLambda = etas2d_spat_intensity(vCat,vX,vY,mU,vPar,fT,varargin)
%
%   Computes the total spatial intensity (Jalilian, 2019; Zhuang et al., 2006)
%   This computes the total intensity directly at each point (x,y) of the grid
%
%   vCat   - event catalogue
%   vX     - vector for X coordinates of the renctanglular grid
%   vY     - vector for Y coordinates of the renctanglular grid
%   mU     - matrix for the background rate u(x,y)
%   vPar   - the parameters of the ETAS model% \mu =vPar(1); A =vPar(2); \alpha =vPar(3); c =vPar(4); p =vPar(5); d =vPar(6); q =vPar(7); \gamma =vPar(8)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 17 May 2021
%   ...
%   version: 1.0.0, 17 May 2021
%
    sPointProc = 'ETAS2D7P'; % 
    for k = 1:length(varargin)
        if strcmp('PointProcess',varargin{k})
            sPointProc = varargin{k+1};
        end
    end

%     d = vPar(6);
%     if strcmp(sPointProc,'ETAS2D5')
%         d = vPar(6)*exp(0.5*vPar(8)*(Parent.mag - fM0));     % for Model 5
%     end
    d2      = vPar(6)^2;
    nX      = length(vX);
    nY      = length(vY);
    mLambda = zeros(nY,nX);
    fF      = 1.0/fT*vPar(2)*(vPar(7)-1)/(pi*d2); % 1/T * A*(q-1)/(\pi*d^2)
    
    nJ = find(vCat(:,1) <= fT,1,'last'); % find the last index for which vCat(:,1) <= vT(j)
    for iy = 1:nY
        for ix = 1:nX
            mLambda(iy,ix) = vPar(1)*mU(iy,ix) + fF*sum( exp(vPar(3).*vCat(1:nJ,4)).*...
                              ((1 + ( (vX(ix) - vCat(1:nJ,3)).^2 + (vY(iy) - vCat(1:nJ,2)).^2 )/d2).^(-vPar(7))) ); % \mu + sum
%             if mLambda(iy,ix) < mU(iy,ix)
%                 disp([mLambda(iy,ix) mU(iy,ix)])
%             end
        end
    end
end



