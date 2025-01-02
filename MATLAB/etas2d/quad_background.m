function fUquadR = quad_background(inR,nJ,vBgProb,vHj,vRk,vRk2,vDk,nNk,varargin)
%
%   Computes the quadrature of the background intensity over a polygon vReg_poly
%
%   vX, vY - linear coordinate vectors over which mU is defined
%   mU     - the function to integrate given as a matrix of size [length(vY) length(vX)]
%   vReg_poly - domain of integration as a polygon
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 11 June 2021
%   ...
%   version: 1.0.0, 15 June 2021
%
    sKernel = 'Gaussian';   % 'Cauchy';   %
    for k = 1:length(varargin)
        if strcmp('Kernel',varargin{k})
            sKernel = varargin{k+1};
        end
    end
    
    f12pi  = 1.0/(2*pi);
    fUquadR = 0;
    for j = 1:nJ
        fIj = 0;
        if inR(j)      % for events inside the target region R
            for n = 1:nNk
                fIj = fIj + vBgProb(j)*f12pi*vDk(j,n)*gauss_kernel_cum(vRk(j,n),vHj(j));
            end
        else
            for n = 1:nNk
                fIR  = vBgProb(j)*f12pi*vDk(j,n)*gauss_kernel_cum(vRk(j,n),vHj(j));
                fIR2 = vBgProb(j)*f12pi*vDk(j,n)*gauss_kernel_cum(vRk2(j,n),vHj(j));
                fIj = fIj + fIR2 - fIR;
            end
        end
        fUquadR = fUquadR + fIj;
    end
end
