function cLegend = plot_gmm_fit(gmm,vXmodel,varargin)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 1.0.1, 29 August 2024
%
    sColor = '#0072BD';
    bComponents = false;
    for k = 1:length(varargin)
        if strcmp('Color',varargin{k})
            sColor = varargin{k+1};
        end
        if strcmp('Components',varargin{k})
            bComponents = true;
        end
    end
    gmm_pdf = pdf(gmm,vXmodel);
    plot(vXmodel,gmm_pdf,'Color',sColor,'LineWidth',1.5);
    hold on
    cLegendGMM = 'GMM: ';
    for n = 1:gmm.NumComponents
        plot(gmm.mu(n),pdf(gmm,gmm.mu(n)),'o','MarkerEdgeColor',sColor,'MarkerFaceColor',sColor,'MarkerSize',6);
        cLegendGMM = [cLegendGMM,'\mu_',num2str(n),' = ',num2str(gmm.mu(n),4)];
        if n < gmm.NumComponents, cLegendGMM = [cLegendGMM,', ']; end
    end
    str = ['AIC = ',num2str(gmm.AIC),', BIC = ',num2str(gmm.BIC)];
    cLegend{1} = sprintf('%s\n%s',cLegendGMM,str);
    cLegend{2} = 'GMM component mean';
    if bComponents && gmm.NumComponents > 1
        for n = 1:gmm.NumComponents
            cLegend{n+1} = '';
        end
        nl = gmm.NumComponents + 1;
        sComColor = {'#4DBEEE','#EDB120','#D95319','#77AC30'};
        for n = 1:gmm.NumComponents
            w = gmm.ComponentProportion(n);
            plot(vXmodel,w*normpdf(vXmodel,gmm.mu(n),sqrt(gmm.Sigma(n))),'Color',sComColor{n},'LineWidth',1.2);
            nl = nl + 1;
            cLegend{nl} = ['Component ',num2str(n),' weight, w = ',num2str(w)];
            plot(gmm.mu(n),pdf(gmm,gmm.mu(n)),'o','MarkerEdgeColor',sComColor{n},'MarkerFaceColor',sComColor{n},'MarkerSize',6);
            nl = nl + 1;
            cLegend{nl} = '';
        end
    end
end

