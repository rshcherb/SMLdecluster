function save_cf(cf,sName,varargin)
%
%   Save current figure into a image file
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 27 April 2023
%   ...
%   version 1.1.0, 4 November 2024
%
    sResolution = '-r600';
    sBestFit = 'off'; % ''
    bAll = true;
    bPdf = false;
    bPng = false;
    bJpeg = false;
    bEps = false;
    bSvg = false;
    bFig = false;
    for k = 1:length(varargin)
        if strcmp('Resolution',varargin{k})
            sResolution = varargin{k+1};
        end
        if strcmp('BestFit',varargin{k})
            sBestFit = varargin{k+1};
        end
        if strcmp('pdf',varargin{k})
            bPdf = true;
            bAll = false;
        end
        if strcmp('png',varargin{k})
            bPng = true;
            bAll = false;
        end
        if strcmp('jpeg',varargin{k})
            bJpeg = true;
            bAll = false;
        end
        if strcmp('eps',varargin{k})
            bEps = true;
            bAll = false;
        end
        if strcmp('svg',varargin{k})
            bSvg = true;
            bAll = false;
        end
        if strcmp('fig',varargin{k})
            bFig = true;
            bAll = false;
        end
    end
    set(cf,'Color','w');
    if bAll || bFig
        savefig(cf,[sName,'.fig'],'compact');
    end
    if bAll || bPng
        print(cf,[sName,'.png'],'-vector','-dpng',sResolution);
    end
    if bAll || bJpeg
        print(cf,[sName,'.jpg'],'-vector','-djpeg',sResolution);
    end
    if bAll || bPdf
        if strcmp(sBestFit,'on')
            print(cf,[sName,'.pdf'],'-vector','-dpdf',sResolution,'-bestfit');
        else
            print(cf,[sName,'.pdf'],'-vector','-dpdf',sResolution);
        end
    end
    if bAll || bEps
        print(cf,[sName,'.eps'],'-vector','-depsc',sResolution);
    end
    if bAll || bSvg
        print(cf,[sName,'.svg'],'-vector','-dsvg',sResolution);
    end
end
