function model_freq_mag(vTM,Model,sTitle,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 17 MArch, 2019
%   ...
%   version 1.1.0, 21 November, 2019
%
    sFitModel   = 'GR'; % 'GR_DTRUNC', 
    fAlpha      = 0.05;
    sPlotStatus = 'on';
    for k = 1:length(varargin)
        if strcmp('FitModel',varargin{k})
            sFitModel = varargin{k+1};
        end
        if strcmp('Alpha',varargin{k})
            fAlpha = varargin{k+1};
        end
        if strcmp('Plot',varargin{k})
            sPlotStatus = varargin{k+1};
        end
    end
    
    if strcmp(sFitModel,'GR')
        % compute GR parameters and plot frequency-magnitude statistics and confidence intervals at fAlpha level
    %     [vParEst, vParErr] = grvalues(vTM(:,2),[Model.fMc, Model.fMmax, Model.fDm],fAlpha);
        [vParEst, vParErr] = gr_fit(vTM(:,2),Model.fMc,Model.fMmax,Model.fDm,'Alpha',fAlpha);
%         [vParEst, vParErr] = gr_fit(vTM(:,2),Model.fMc,Model.fMmax,Model.fDm,'Alpha',fAlpha,...
%             'FitMethod','Utsu','ErrorMethod','Aki','Display','on');
%             'FitMethod','Utsu','ErrorMethod','ShiBolt','Display','on');
%             'FitMethod','Bender','ErrorMethod','ShiBolt','Display','on');
%             'FitMethod','Bender','ErrorMethod','Aki','Display','on');
    elseif strcmp(sFitModel,'GR_DTRUNC')
    end
    
    if strcmp(sPlotStatus,'on')
        figure('Name','GR scaling');
        %
        set(gcf,'color','w'); % background color for the figure
        % plot frequency-magnitude distribution and GR scaling with given b
        [vXbin, vBinData, vCumData] = gr_fit_plot(vTM(:,2),vParEst,vParErr,Model.fMc,Model.fDm,'Alpha',fAlpha,...
                                                  'PoissonBounds','BathLaw');
        title(sTitle);
        hold off;
        
        sName = [Model.sFileName,'_gr'];
        mTmp = [vXbin', vBinData', vCumData'];
        save([sName,'.dat'],'mTmp','-ascii');
        saveas(gcf,[sName,'.png']);
        saveas(gcf,[sName,'.pdf']);
    end
end

