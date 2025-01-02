function mFBs = fbsurf(nMax,fSigma,fH,bAdd)
%
% Based on the midpoin displacement and successive additions algorithm from
% Peitgen H.-O., Saupe D. (eds.) The Science of Fractal Images, Springer, 1988, p. 100.
%
% Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%   version 1.0.0, November 26, 2014
%
    N = 2^nMax;
    delta = fSigma;
    mFBs = zeros(N+1,N+1);
    mFBs(1,1)     = delta*randn();
    mFBs(1,N+1)   = delta*randn();
    mFBs(N+1,1)   = delta*randn();
    mFBs(N+1,N+1) = delta*randn();
    
    nD = N; 
    nd = N/2;    
    for n = 1:nMax
        
        % grid type I to type II
        delta = delta*0.5^(0.5*fH);
        
        x = nd+1:nD:N+1-nd;
        lenx = length(x);
        mFBs(x,x) = (mFBs(x+nd,x+nd) + mFBs(x-nd,x+nd) + mFBs(x+nd,x-nd) + mFBs(x-nd,x-nd))/4.0 +...
            delta*randn(lenx,lenx);

        % displace other points
        if bAdd
            x1 = 1:nD:N+1;
            lenx1 = length(x1);
            mFBs(x1,x1) = mFBs(x1,x1) + delta*randn(lenx1);
        end
        
        % grid type II to type I
        delta = delta*0.5^(0.5*fH);
        
        % boundary grid points
        mFBs(x,1)   = (mFBs(x+nd,1)   + mFBs(x-nd,1)   + mFBs(x,nd))/3.0     + delta*randn(lenx,1);
        mFBs(x,N+1) = (mFBs(x+nd,N+1) + mFBs(x-nd,N+1) + mFBs(x,N+1-nd))/3.0 + delta*randn(lenx,1);
        mFBs(1,x)   = (mFBs(1,x+nd)   + mFBs(1,x-nd)   + mFBs(nd,x))/3.0     + delta*randn(1,lenx);
        mFBs(N+1,x) = (mFBs(N+1,x+nd) + mFBs(N+1,x-nd) + mFBs(N+1-nd,x))/3.0 + delta*randn(1,lenx);

        % interior grid points
        y = nD+1:nD:N+1-nd;
        leny = length(y);
        mFBs(x,y) = (mFBs(x+nd,y) + mFBs(x-nd,y) + mFBs(x,y+nd) + mFBs(x,y-nd))/4.0 + delta*randn(lenx,leny);
        mFBs(y,x) = (mFBs(y+nd,x) + mFBs(y-nd,x) + mFBs(y,x+nd) + mFBs(y,x-nd))/4.0 + delta*randn(leny,lenx);
        
        % displace other points
        if bAdd
            mFBs(x,x) = mFBs(x,x) + delta*randn(lenx);
            mFBs(y,y) = mFBs(y,y) + delta*randn(leny);
        end
        nd = nd/2;
        nD = nD/2;
    end
end

